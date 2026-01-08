// -----------------------------------------------------------------------
// Copyright 2018 Autodesk, Inc. All rights reserved.
// 
// Use of this software is subject to the terms of the Autodesk license
// agreement provided at the time of installation or download, or which
// otherwise accompanies this software in either electronic or hard copy form.
// -----------------------------------------------------------------------

using FeatureCAM;
using Microsoft.Win32;
using Newtonsoft.Json;
using System;
using System.Collections.Generic;
using System.Data;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Runtime.InteropServices;
using System.Runtime.InteropServices.ComTypes;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;


namespace FeatureCAMToCAMplete
{

    public class FCToCAMplete
    {
        private static Queue<string> _batchQueue;
        private static string _batchOutDir;
        private static Timer _batchTimer;
        static public UI main_form = null;

        public class CncDfMConfig
        {
            public double coplanarEpsDeg = 1.0;
            public double sharpEdgeAngleDeg = 80.0;

            // Units should match your STL export units (your bbox suggests inches)
            public double minEndmillDiameter = 0.0625;   // 1/16"
            public double defaultEndmillDiameter = 0.25; // fallback if radius unknown

            // STL sampling resolution for ray casting (trade speed vs accuracy)
            public int rayGridNx = 70;
            public int rayGridNy = 70;

            // Classifiers / thresholds
            public double sharpEdgeDihedralThresholdDeg = 120.0; // edge is "sharp" if dihedral < this
            public int sharpEdgeCountForSharpCorners = 50;
            public double minDihedralForSharpCornersDeg = 95.0;

            // Scoring
            public double depthToDiamGood = 4.0;
            public double depthToDiamWarn = 8.0;
            public double minTopAccessibilityPctGood = 0.90; // 90%
            public double minTopAccessibilityPctWarn = 0.75; // 75%
        }

        public class CncMetrics
        {
            public string partName { get; set; }

            public int solids { get; set; }
            public int setupsInFile { get; set; }     // doc.Setups.Count (not "estimated")
            public int operations { get; set; }
            public int tools { get; set; }

            // Dimensions from STL bbox
            public double partX { get; set; }
            public double partY { get; set; }
            public double partZ { get; set; }

            public double stockX { get; set; }
            public double stockY { get; set; }
            public double stockZ { get; set; }

            // Mesh/DFM metrics
            public double estimatedPocketDepth { get; set; }       // from STL ray casting
            public double estimatedMinToolDiameter { get; set; }   // from sharp corner proxy
            public double depthToToolDiamRatio { get; set; }       // pocketDepth / minToolDiam

            public double accessibleFromTopPct { get; set; }       // visibility of upward-facing surfaces from +Z
            public int estimatedSetupsMin { get; set; }            // 1 = top only, 2 = top + flip

            // Sharp internal corner proxy
            public double? minInternalRadius { get; set; }         // null unknown, 0 sharp
            public int sharedEdgeCount { get; set; }
            public int sharpEdgeCount { get; set; }
            public double minDihedralDeg { get; set; }
            public double concavitySharpnessScore { get; set; }
            // Backwards-compat aliases (for older code that expects these names)
            // Backwards-compat aliases (old code can still assign)
            [JsonIgnore]
            public int setups
            {
                get => estimatedSetupsMin;
                set => estimatedSetupsMin = value;
            }

            [JsonIgnore]
            public double maxPocketDepth
            {
                get => estimatedPocketDepth;
                set => estimatedPocketDepth = value;
            }

            [JsonIgnore]
            public double depthToWidthRatio
            {
                get => depthToToolDiamRatio;
                set => depthToToolDiamRatio = value;
            }



            // Output
            public double viabilityScore { get; set; }
            public List<string> flags { get; set; } = new List<string>();

            public string ToJson() => JsonConvert.SerializeObject(this, Formatting.Indented);
        }
        private static bool IsFinite(double x)
        {
            return !double.IsNaN(x) && !double.IsInfinity(x);
        }

        public static CncMetrics ComputeCncMetrics(FeatureCAM.FMDocument doc, string outputDir, CncDfMConfig cfg = null)
        {
            cfg = cfg ?? new CncDfMConfig();

            // Clamp config so we never get 0 tool diam
            cfg.minEndmillDiameter = Math.Max(cfg.minEndmillDiameter, 1e-6);
            cfg.defaultEndmillDiameter = Math.Max(cfg.defaultEndmillDiameter, cfg.minEndmillDiameter);

            var m = new CncMetrics();
            m.flags = m.flags ?? new List<string>();

            // Basic doc stats
            try { m.partName = doc.PartName ?? ""; } catch { m.partName = ""; }
            try { m.solids = doc.Solids != null ? doc.Solids.Count : 0; } catch { m.solids = 0; }
            try { m.setupsInFile = doc.Setups != null ? doc.Setups.Count : 0; } catch { m.setupsInFile = 0; }
            try { m.operations = doc.Operations != null ? doc.Operations.Count : 0; } catch { m.operations = 0; }

            // If you have CountUniqueTools, keep it; else leave tools=0
            try { m.tools = CountUniqueTools(doc); } catch { m.tools = 0; }

            Directory.CreateDirectory(outputDir);

            // --- Export + analyze PART STL ---
            if (m.solids < 1)
            {
                m.flags.Add("No solids in document");
                m.viabilityScore = 0;
                return m;
            }

            string partStl = Path.Combine(outputDir, "_tmp_part.stl");

            try
            {
                var solid1 = (FeatureCAM.FMSolid)doc.Solids.Item(1);
                if (!ExportToStl(solid1, partStl, out string err))
                {
                    m.flags.Add("ExportToSTL(part) failed: " + err);
                    m.viabilityScore = 0;
                    return m;
                }
            }
            catch (Exception ex)
            {
                m.flags.Add("ExportToSTL(part) exception: " + ex.Message);
                m.viabilityScore = 0;
                return m;
            }

            // Part bbox from STL (ASCII STL expected)
            if (TryComputeBboxFromAsciiStl(partStl, out double pminx, out double pminy, out double pminz,
                                                     out double pmaxx, out double pmaxy, out double pmaxz))
            {
                m.partX = pmaxx - pminx;
                m.partY = pmaxy - pminy;
                m.partZ = pmaxz - pminz;
                DebugProbe.Touch($"STL bbox part => [{pminx},{pminy},{pminz}]..[{pmaxx},{pmaxy},{pmaxz}]");
            }
            else
            {
                m.flags.Add("Part STL bbox parse failed (binary STL or unexpected format)");
            }

            // Concavity / sharp edges (dihedral) from STL
            try
            {
                AnalyzeSharpEdgesFromStl(
                    partStl,
                    cfg,
                    out int sharedEdges,
                    out int sharpEdges,
                    out double minDihedralDeg,
                    out double sharpScore);

                m.sharedEdgeCount = sharedEdges;
                m.sharpEdgeCount = sharpEdges;
                m.minDihedralDeg = minDihedralDeg;
                m.concavitySharpnessScore = sharpScore;

                DebugProbe.Touch($"Concavity: sharedEdges={sharedEdges} sharpEdges={sharpEdges} minDihedral={minDihedralDeg:F1} sharpScore={sharpScore:F4}");

                // Interpret sharp internal corners as radius ~ 0 (proxy)
                // (This is heuristic; later we can estimate a real radius from mesh curvature)
                if (sharpEdges > cfg.sharpEdgeCountForSharpCorners && minDihedralDeg < cfg.minDihedralForSharpCornersDeg)
                    m.minInternalRadius = 0.0;
                else
                    m.minInternalRadius = null;
            }
            catch (Exception ex)
            {
                m.flags.Add("Concavity analysis failed: " + ex.Message);
            }

            // Estimated min tool diameter: if sharp internal corners, tool must be <= shop min (proxy)
            if (m.minInternalRadius.HasValue && m.minInternalRadius.Value <= 1e-9)
                m.estimatedMinToolDiameter = cfg.minEndmillDiameter;
            else
                m.estimatedMinToolDiameter = cfg.defaultEndmillDiameter;

            // Clamp & log
            m.estimatedMinToolDiameter = Math.Max(m.estimatedMinToolDiameter, 1e-6);
            DebugProbe.Touch($"MinToolDiam={m.estimatedMinToolDiameter} (cfg.min={cfg.minEndmillDiameter}, cfg.def={cfg.defaultEndmillDiameter})");

            // Ray/mesh metrics (accessibility + pocket depth)
            try
            {
                // IMPORTANT: this must be the version where normals are computed from vertices
                // (ReadAsciiStlTris computes normals via cross product)
                var tris = ReadAsciiStlTris(partStl);

                ComputeTopAccessibilityAndPocketDepth(tris, cfg,
                    out double accessPct,
                    out double pocketDepth);

                m.accessibleFromTopPct = accessPct;
                m.estimatedPocketDepth = pocketDepth;

                DebugProbe.Touch($"TopAccess={m.accessibleFromTopPct:F3} PocketDepth={m.estimatedPocketDepth}");
            }
            catch (Exception ex)
            {
                m.flags.Add("Top accessibility / pocket depth failed: " + ex.Message);
            }

            // Depth/tool ratio (main driver)
            if (m.estimatedMinToolDiameter > 1e-9)
                m.depthToToolDiamRatio = m.estimatedPocketDepth / m.estimatedMinToolDiameter;
            else
                m.depthToToolDiamRatio = 0.0;

            // Estimated setups (3-axis): heuristic
            // If there is meaningful downward-facing geometry => likely need flip
            try
            {
                var tris2 = ReadAsciiStlTris(partStl);
                m.estimatedSetupsMin = HasMeaningfulDownFacing(tris2) ? 2 : 1;
            }
            catch
            {
                m.estimatedSetupsMin = 1;
            }

            // --- Export + analyze STOCK STL (optional) ---
            try
            {
                object stockObj = null;
                try { stockObj = doc.Stock; } catch { stockObj = null; }

                if (stockObj != null)
                {
                    string stockStl = Path.Combine(outputDir, "_tmp_stock.stl");
                    if (ExportToStl(stockObj, stockStl, out string errStock))
                    {
                        if (TryComputeBboxFromAsciiStl(stockStl, out double sminx, out double sminy, out double sminz,
                                                                 out double smaxx, out double smaxy, out double smaxz))
                        {
                            m.stockX = smaxx - sminx;
                            m.stockY = smaxy - sminy;
                            m.stockZ = smaxz - sminz;
                            DebugProbe.Touch($"STL bbox stock => [{sminx},{sminy},{sminz}]..[{smaxx},{smaxy},{smaxz}]");
                        }
                    }
                }
            }
            catch { /* stock is optional */ }

            // Score + flags
            try
            {
                m.viabilityScore = Score3Axis(m, cfg);
            }
            catch
            {
                m.viabilityScore = 0;
                m.flags.Add("Scoring failed");
            }
            bool impossibleSharpInsideCorners = (m.minInternalRadius.HasValue && m.minInternalRadius.Value <= 1e-9);

            if (impossibleSharpInsideCorners)
            {
                m.flags.Add("IMPOSSIBLE (3-axis endmill): sharp internal corners without relief/EDM");
            }

            return m;
        }


        // ----------------------------
        // 3-axis score (tool ratio, accessibility, setups)
        // ----------------------------
        private static double Score3Axis(CncMetrics m, CncDfMConfig cfg)
        {
            double score = 100.0;

            // Tool reach/stickout proxy (most important)
            if (m.depthToToolDiamRatio > cfg.depthToDiamGood)
            {
                score -= 10;
                m.flags.Add($"Depth/ToolDiam > {cfg.depthToDiamGood:F1}");
            }
            if (m.depthToToolDiamRatio > cfg.depthToDiamWarn)
            {
                score -= 25;
                m.flags.Add($"Depth/ToolDiam > {cfg.depthToDiamWarn:F1} (reach risk)");
            }

            // Sharp internal corners (forces smaller tools / more ops)
            if (m.minInternalRadius.HasValue && m.minInternalRadius.Value <= 1e-9)
            {
                score -= 15;
                m.flags.Add("Sharp internal corners (small tool / extra ops)");
            }

            // Accessibility from top (3-axis)
            if (m.accessibleFromTopPct < cfg.minTopAccessibilityPctGood)
            {
                score -= 10;
                m.flags.Add("Top accessibility < 90%");
            }
            if (m.accessibleFromTopPct < cfg.minTopAccessibilityPctWarn)
            {
                score -= 20;
                m.flags.Add("Top accessibility < 75% (risk of non-3axis features)");
            }

            // Setups (fixturing/time)
            if (m.estimatedSetupsMin >= 2)
            {
                score -= 10;
                m.flags.Add("Likely requires flip (>=2 setups)");
            }
            if (m.flags.Contains("IMPOSSIBLE (3-axis endmill): sharp internal corners without relief/EDM"))
            {
                // hard cap; or set to 0 if you prefer
                score = Math.Min(score, 20);
            }

            // Clamp
            if (score < 0) score = 0;
            if (score > 100) score = 100;
            if (m.minInternalRadius.HasValue && m.minInternalRadius.Value <= 1e-9)
            {
                m.flags.Add("IMPOSSIBLE (3-axis): sharp internal corners (needs relief or EDM)");
                score = Math.Min(score, 20); // or score = 0;
            }

            return score;
        }

        // ============================================================
        // Mesh analysis: pocket depth + top accessibility (3-axis +Z)
        // ============================================================

        // Simple vector/triangle structs
        private struct Vec3
        {
            public double X, Y, Z;
            public Vec3(double x, double y, double z) { X = x; Y = y; Z = z; }
            public static Vec3 operator -(Vec3 a, Vec3 b) => new Vec3(a.X - b.X, a.Y - b.Y, a.Z - b.Z);
            public static double Dot(Vec3 a, Vec3 b) => a.X * b.X + a.Y * b.Y + a.Z * b.Z;
            public static Vec3 Cross(Vec3 a, Vec3 b) =>
                new Vec3(a.Y * b.Z - a.Z * b.Y, a.Z * b.X - a.X * b.Z, a.X * b.Y - a.Y * b.X);
            public double Norm() => Math.Sqrt(X * X + Y * Y + Z * Z);
            public Vec3 Normalize()
            {
                var n = Norm();
                return (n < 1e-12) ? new Vec3(0, 0, 0) : new Vec3(X / n, Y / n, Z / n);
            }
        }

        private struct Tri
        {
            public Vec3 N;
            public Vec3 A, B, C;
        }

        private static List<Tri> ReadAsciiStlTris(string stlPath)
        {
            var tris = new List<Tri>();
            using (var sr = new StreamReader(stlPath))
            {
                string line;
                Vec3 a = new Vec3(), b = new Vec3(), c = new Vec3();
                int vcount = 0;

                while ((line = sr.ReadLine()) != null)
                {
                    var s = line.Trim();

                    if (s.StartsWith("vertex", StringComparison.OrdinalIgnoreCase))
                    {
                        var parts = s.Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                        if (parts.Length >= 4)
                        {
                            var v = new Vec3(
                                double.Parse(parts[1], CultureInfo.InvariantCulture),
                                double.Parse(parts[2], CultureInfo.InvariantCulture),
                                double.Parse(parts[3], CultureInfo.InvariantCulture));

                            if (vcount == 0) a = v;
                            else if (vcount == 1) b = v;
                            else if (vcount == 2) c = v;

                            vcount++;
                        }
                    }
                    else if (s.StartsWith("endfacet", StringComparison.OrdinalIgnoreCase))
                    {
                        if (vcount >= 3)
                        {
                            // Compute normal from vertices
                            var n = Vec3.Cross(b - a, c - a).Normalize();
                            tris.Add(new Tri { N = n, A = a, B = b, C = c });
                        }
                        vcount = 0;
                    }
                }
            }
            return tris;
        }



        // Möller–Trumbore ray/triangle intersection
        // Ray: origin O, direction D (assume D = (0,0,-1))
        private static bool RayIntersectTri(Vec3 O, Vec3 D, Tri t, out double hitT)
        {
            hitT = 0;
            const double EPS = 1e-9;

            var e1 = t.B - t.A;
            var e2 = t.C - t.A;
            var p = Vec3.Cross(D, e2);
            double det = Vec3.Dot(e1, p);
            if (Math.Abs(det) < EPS) return false;
            double invDet = 1.0 / det;

            var s = O - t.A;
            double u = invDet * Vec3.Dot(s, p);
            if (u < -EPS || u > 1.0 + EPS) return false;

            var q = Vec3.Cross(s, e1);
            double v = invDet * Vec3.Dot(D, q);
            if (v < -EPS || u + v > 1.0 + EPS) return false;

            double tt = invDet * Vec3.Dot(e2, q);
            if (tt > EPS)
            {
                hitT = tt;
                return true;
            }
            return false;
        }

        private static void ComputeTopAccessibilityAndPocketDepth(
    List<Tri> tris,
    CncDfMConfig cfg,
    out double accessibleFromTopPct,
    out double estimatedPocketDepth)
        {
            cfg.rayGridNx = 250;
            cfg.rayGridNy = 250;

            accessibleFromTopPct = 0.0;
            estimatedPocketDepth = 0.0;

            if (tris == null || tris.Count == 0) return;

            // bbox
            double minx = double.PositiveInfinity, miny = double.PositiveInfinity, minz = double.PositiveInfinity;
            double maxx = double.NegativeInfinity, maxy = double.NegativeInfinity, maxz = double.NegativeInfinity;
            for (int i = 0; i < tris.Count; i++)
            {
                Update(ref minx, ref miny, ref minz, ref maxx, ref maxy, ref maxz, tris[i].A);
                Update(ref minx, ref miny, ref minz, ref maxx, ref maxy, ref maxz, tris[i].B);
                Update(ref minx, ref miny, ref minz, ref maxx, ref maxy, ref maxz, tris[i].C);
            }

            double dx = maxx - minx, dy = maxy - miny, dz = maxz - minz;
            if (dx < 1e-9 || dy < 1e-9) return;

            int nx = Math.Max(30, cfg.rayGridNx);
            int ny = Math.Max(30, cfg.rayGridNy);

            var D = new Vec3(0, 0, -1);
            double zStart = maxz + 5.0 * Math.Max(dx, dy) + 1.0;

            // eps for intersection distinctness
            double diag = Math.Sqrt(dx * dx + dy * dy + dz * dz);
            double epsZ = Math.Max(1e-4, 1e-6 * diag);

            // Collect first-hit Z for all rays that hit something
            var firstHits = new List<double>(nx * ny);

            int total = nx * ny;
            int hitAny = 0;

            // temp list for hits along one ray
            var hits = new List<double>(16);

            for (int iy = 0; iy < ny; iy++)
            {
                double y = miny + (iy + 0.5) * dy / ny;
                for (int ix = 0; ix < nx; ix++)
                {
                    double x = minx + (ix + 0.5) * dx / nx;
                    var O = new Vec3(x, y, zStart);

                    hits.Clear();

                    for (int t = 0; t < tris.Count; t++)
                    {
                        if (RayIntersectTri(O, D, tris[t], out double tt))
                        {
                            double zHit = zStart - tt;
                            hits.Add(zHit);
                        }
                    }

                    if (hits.Count == 0)
                        continue;

                    hitAny++;

                    // "first hit from above" = max zHit
                    double zFirst = hits[0];
                    for (int k = 1; k < hits.Count; k++)
                        if (hits[k] > zFirst) zFirst = hits[k];

                    firstHits.Add(zFirst);
                }
            }

            accessibleFromTopPct = (total > 0) ? ((double)hitAny / total) : 0.0;

            if (firstHits.Count == 0)
            {
                estimatedPocketDepth = 0.0;
                DebugProbe.Touch($"TopAccess raysHit={hitAny}/{total} => {accessibleFromTopPct:F3}, pocketDepth=0 (no hits)");
                return;
            }

            // Sort first-hit heights to estimate top surface height robustly
            firstHits.Sort(); // ascending
                              // Use a high percentile (e.g. 98th) to avoid a single spike / noise
            int idxTop = (int)Math.Round(0.98 * (firstHits.Count - 1));
            if (idxTop < 0) idxTop = 0;
            if (idxTop >= firstHits.Count) idxTop = firstHits.Count - 1;

            double zTop = firstHits[idxTop];

            // Identify "below-top" first hits as candidate pocket floors
            // Threshold: must be meaningfully below zTop
            // Use 2% of thickness OR epsZ, whichever is larger.
            double thresh = Math.Max(epsZ * 10.0, 0.02 * dz);

            double bestPocketFloorZ = double.NegativeInfinity;
            int pocketRayCount = 0;

            for (int i = 0; i < firstHits.Count; i++)
            {
                double z = firstHits[i];
                if (zTop - z > thresh)
                {
                    pocketRayCount++;
                    if (z > bestPocketFloorZ) bestPocketFloorZ = z;
                }
            }

            if (pocketRayCount > 0 && bestPocketFloorZ > double.NegativeInfinity / 2)
                estimatedPocketDepth = zTop - bestPocketFloorZ;
            else
                estimatedPocketDepth = 0.0;

            DebugProbe.Touch(
                $"TopAccess raysHit={hitAny}/{total} => {accessibleFromTopPct:F3}, " +
                $"zTop={zTop}, pocketRays={pocketRayCount}, pocketDepth={estimatedPocketDepth}");
            double pocketMmAssumingInches = estimatedPocketDepth * 25.4;
            DebugProbe.Touch($"PocketDepth(mm assuming inches) = {pocketMmAssumingInches:F2}");

        }





        private static void Update(ref double minx, ref double miny, ref double minz,
                                   ref double maxx, ref double maxy, ref double maxz, Vec3 v)
        {
            if (v.X < minx) minx = v.X; if (v.Y < miny) miny = v.Y; if (v.Z < minz) minz = v.Z;
            if (v.X > maxx) maxx = v.X; if (v.Y > maxy) maxy = v.Y; if (v.Z > maxz) maxz = v.Z;
        }

        private static bool HasMeaningfulDownFacing(List<Tri> tris)
        {
            // If there is significant downward-facing area, you'll likely need a flip to machine bottom features.
            // This is a heuristic. Tune if you want.
            int down = 0;
            for (int i = 0; i < tris.Count; i++)
                if (tris[i].N.Z < -0.2) down++;
            return down > 100; // threshold
        }

        // ============================================================
        // Sharp internal corners proxy (edge dihedral)
        // ============================================================
        private static void AnalyzeSharpEdgesFromStl(
    string stlPath,
    CncDfMConfig cfg,
    out int sharedEdgeCount,
    out int sharpEdgeCount,
    out double minDihedralDeg,
    out double sharpScore)
        {
            sharedEdgeCount = 0;
            sharpEdgeCount = 0;
            minDihedralDeg = 0.0;
            sharpScore = 0.0;

            // Local constants so you DON'T need new cfg fields
            const double COPLANAR_EPS_DEG = 1.0;  // ignore 0–1° coplanar triangulation seams
            const double SHARP_EDGE_DEG = 80.0; // classify >=80° as "sharp"

            // IMPORTANT: this must be the STL reader that recomputes normals from vertices
            List<Tri> tris = ReadAsciiStlTris(stlPath);
            if (tris == null || tris.Count == 0) return;

            // Build edge map: undirected edge -> triangles using it
            var edgeMap = new Dictionary<string, List<int>>(tris.Count * 3);

            for (int i = 0; i < tris.Count; i++)
            {
                AddEdge(edgeMap, i, tris[i].A, tris[i].B);
                AddEdge(edgeMap, i, tris[i].B, tris[i].C);
                AddEdge(edgeMap, i, tris[i].C, tris[i].A);
            }

            int nonCoplanarCount = 0;
            double minNonCoplanar = 180.0;

            foreach (var kv in edgeMap)
            {
                var triIds = kv.Value;
                if (triIds == null || triIds.Count != 2)
                    continue;

                Tri t1 = tris[triIds[0]];
                Tri t2 = tris[triIds[1]];

                Vec3 n1 = t1.N;
                Vec3 n2 = t2.N;

                double dot = Vec3.Dot(n1, n2);
                if (dot > 1.0) dot = 1.0;
                else if (dot < -1.0) dot = -1.0;

                double dihedral = Math.Acos(dot) * 180.0 / Math.PI;

                // Debug first few edges
                if (sharedEdgeCount < 5)
                {
                    DebugProbe.Touch(
                        $"DBG dihedral={dihedral:F2} dot={dot:F4} " +
                        $"n1=({n1.X:F3},{n1.Y:F3},{n1.Z:F3}) " +
                        $"n2=({n2.X:F3},{n2.Y:F3},{n2.Z:F3})");
                }

                sharedEdgeCount++;

                // Ignore coplanar seams (these force min=0 otherwise)
                if (dihedral <= COPLANAR_EPS_DEG)
                    continue;

                nonCoplanarCount++;

                if (dihedral < minNonCoplanar)
                    minNonCoplanar = dihedral;

                if (dihedral >= SHARP_EDGE_DEG)
                    sharpEdgeCount++;
            }

            minDihedralDeg = (nonCoplanarCount > 0) ? minNonCoplanar : 0.0;
            sharpScore = (nonCoplanarCount > 0) ? ((double)sharpEdgeCount / nonCoplanarCount) : 0.0;
        }

        // --- Helpers (put these in the SAME class as AnalyzeSharpEdgesFromStl) ---

        private static void AddEdge(Dictionary<string, List<int>> edgeMap, int triIndex, Vec3 v1, Vec3 v2)
        {
            // Quantize to match vertices robustly
            const double Q = 1e-5;

            string k1 = VKey(v1, Q);
            string k2 = VKey(v2, Q);

            // undirected => sort endpoints
            string key = (string.CompareOrdinal(k1, k2) <= 0) ? (k1 + "__" + k2) : (k2 + "__" + k1);

            List<int> list;
            if (!edgeMap.TryGetValue(key, out list))
            {
                list = new List<int>(2);
                edgeMap[key] = list;
            }
            list.Add(triIndex);
        }

        private static string VKey(Vec3 v, double q)
        {
            long xi = (long)Math.Round(v.X / q);
            long yi = (long)Math.Round(v.Y / q);
            long zi = (long)Math.Round(v.Z / q);
            return xi.ToString() + "|" + yi.ToString() + "|" + zi.ToString();
        }

        // ============================================================
        // FeatureCAM helpers (tools + STL export + bbox parsing)
        // ============================================================
        private static int CountUniqueTools(FeatureCAM.FMDocument doc)
        {
            var keys = new HashSet<string>(StringComparer.OrdinalIgnoreCase);
            int opCount = SafeInt(() => doc.Operations.Count);

            for (int i = 1; i <= opCount; i++)
            {
                try
                {
                    var op = (FeatureCAM.FMOperation)doc.Operations.Item(i);
                    object toolObj = null;
                    try { toolObj = op.Tool; } catch { }
                    string key = GetToolKey(toolObj);
                    if (!string.IsNullOrEmpty(key)) keys.Add(key);
                }
                catch { }
            }

            return keys.Count;
        }

        private static string GetToolKey(object toolObj)
        {
            if (toolObj == null) return "";
            foreach (var p in new[] { "ToolNumber", "Number", "Name", "Description" })
            {
                try
                {
                    object v = toolObj.GetType().InvokeMember(p, BindingFlags.GetProperty, null, toolObj, null);
                    if (v != null)
                    {
                        var s = v.ToString();
                        if (!string.IsNullOrWhiteSpace(s)) return p + ":" + s;
                    }
                }
                catch { }
            }
            return "Type:" + toolObj.GetType().FullName;
        }

        private static bool ExportToStl(object comObj, string stlPath, out string err)
        {
            err = "";
            try
            {
                // ExportToSTL(string path, out string err)
                object[] args = new object[] { stlPath, "" };
                comObj.GetType().InvokeMember("ExportToSTL", BindingFlags.InvokeMethod, null, comObj, args);
                err = (args.Length > 1 && args[1] != null) ? args[1].ToString() : "";
                return string.IsNullOrEmpty(err);
            }
            catch (Exception ex)
            {
                err = ex.Message;
                return false;
            }
        }

        private static bool TryComputeBboxFromAsciiStl(
            string stlPath,
            out double minx, out double miny, out double minz,
            out double maxx, out double maxy, out double maxz)
        {
            minx = miny = minz = double.PositiveInfinity;
            maxx = maxy = maxz = double.NegativeInfinity;
            if (!File.Exists(stlPath)) return false;

            bool any = false;
            foreach (var line in File.ReadLines(stlPath))
            {
                var s = line.Trim();
                if (!s.StartsWith("vertex ", StringComparison.OrdinalIgnoreCase)) continue;

                var parts = s.Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                if (parts.Length < 4) continue;

                if (double.TryParse(parts[1], NumberStyles.Float, CultureInfo.InvariantCulture, out var x) &&
                    double.TryParse(parts[2], NumberStyles.Float, CultureInfo.InvariantCulture, out var y) &&
                    double.TryParse(parts[3], NumberStyles.Float, CultureInfo.InvariantCulture, out var z))
                {
                    any = true;
                    if (x < minx) minx = x; if (y < miny) miny = y; if (z < minz) minz = z;
                    if (x > maxx) maxx = x; if (y > maxy) maxy = y; if (z > maxz) maxz = z;
                }
            }
            return any && !double.IsInfinity(minx);
        }

        private static int SafeInt(Func<int> f) { try { return f(); } catch { return 0; } }

        static public FeatureCAM.Application Application
        {
            get { return fc; }
            set { fc = value; }
        }
        static private FeatureCAM.Application fc;

        public FCToCAMplete() { }
        private static void LogDocCreationSignatures()
        {
            object docsObj = fc.Documents;

            var fm = GetComParamVarTypes(docsObj, "AddFM");
            DebugProbe.Touch("AddFM VTs: " + (fm == null ? "null" : string.Join(",", fm)));

            var mf = GetComParamVarTypes(docsObj, "AddMF");
            DebugProbe.Touch("AddMF VTs: " + (mf == null ? "null" : string.Join(",", mf)));

            var tsf = GetComParamVarTypes(docsObj, "AddTSF");
            DebugProbe.Touch("AddTSF VTs: " + (tsf == null ? "null" : string.Join(",", tsf)));
        }

        private static short[] GetComParamVarTypes(object comObj, string methodName)
        {
            var disp = comObj as IDispatch;
            if (disp == null) return null;

            Guid riid = Guid.Empty;
            int[] dispIds = new int[1];
            disp.GetIDsOfNames(ref riid, new[] { methodName }, 1, 0, dispIds);
            int targetMemId = dispIds[0];

            disp.GetTypeInfo(0, 0, out IntPtr pTypeInfo);
            var typeInfo = (System.Runtime.InteropServices.ComTypes.ITypeInfo)Marshal.GetObjectForIUnknown(pTypeInfo);

            typeInfo.GetTypeAttr(out IntPtr pAttr);
            var attr = (System.Runtime.InteropServices.ComTypes.TYPEATTR)Marshal.PtrToStructure(
                pAttr, typeof(System.Runtime.InteropServices.ComTypes.TYPEATTR));

            try
            {
                for (int i = 0; i < attr.cFuncs; i++)
                {
                    typeInfo.GetFuncDesc(i, out IntPtr pFuncDesc);
                    var funcDesc = (System.Runtime.InteropServices.ComTypes.FUNCDESC)Marshal.PtrToStructure(
                        pFuncDesc, typeof(System.Runtime.InteropServices.ComTypes.FUNCDESC));

                    try
                    {
                        if (funcDesc.memid != targetMemId) continue;

                        int n = funcDesc.cParams;
                        short[] vts = new short[n];

                        // FUNCDESC.lprgelemdescParam points to an array of ELEMDESC
                        int elemSize = Marshal.SizeOf(typeof(System.Runtime.InteropServices.ComTypes.ELEMDESC));
                        for (int p = 0; p < n; p++)
                        {
                            IntPtr pElem = new IntPtr(funcDesc.lprgelemdescParam.ToInt64() + p * elemSize);
                            var elem = (System.Runtime.InteropServices.ComTypes.ELEMDESC)Marshal.PtrToStructure(
                                pElem, typeof(System.Runtime.InteropServices.ComTypes.ELEMDESC));

                            vts[p] = elem.tdesc.vt; // VARTYPE (VT_BSTR, VT_I4, VT_BOOL, etc.)
                        }
                        return vts;
                    }
                    finally
                    {
                        typeInfo.ReleaseFuncDesc(pFuncDesc);
                    }
                }
                return null;
            }
            finally
            {
                typeInfo.ReleaseTypeAttr(pAttr);
            }
        }
        
        [ComImport, Guid("00020400-0000-0000-C000-000000000046"),
 InterfaceType(ComInterfaceType.InterfaceIsIUnknown)]
        interface IDispatch
        {
            int GetTypeInfoCount();
            void GetTypeInfo(int iTInfo, int lcid, out IntPtr ppTInfo);
            void GetIDsOfNames(ref Guid riid,
                [MarshalAs(UnmanagedType.LPArray, ArraySubType = UnmanagedType.LPWStr)] string[] rgsNames,
                int cNames, int lcid,
                [MarshalAs(UnmanagedType.LPArray)] int[] rgDispId);
            void Invoke(int dispIdMember, ref Guid riid, int lcid, short wFlags,
                IntPtr pDispParams, out object pVarResult,
                IntPtr pExcepInfo, IntPtr puArgErr);
        }

        private static FeatureCAM.Application GetFeatureCamApp(object obj)
        {
            // 1) Try direct cast first
            try
            {
                var app = obj as FeatureCAM.Application;
                if (app != null) return app;
            }
            catch { /* ignore */ }

            // 2) Try COM running object table
            try
            {
                object o = Marshal.GetActiveObject("FeatureCAM.Application");
                var app = o as FeatureCAM.Application;
                if (app != null) return app;
            }
            catch { /* ignore */ }

            // 3) As a last resort, try ProgID variants sometimes used
            try
            {
                object o = Marshal.GetActiveObject("FeatureCAM.Application.1");
                var app = o as FeatureCAM.Application;
                if (app != null) return app;
            }
            catch { /* ignore */ }

            return null;
        }
        public static void RunMetricsTest()
        {
            DebugProbe.Touch("RunMetricsTest invoked");
            MessageBox.Show("RunMetricsTest fired!", "FeatureCAM Add-in");
        }
        public static void RunBatch(string inputDir, string outputDir)
        {
            LogDocCreationSignatures();
            DebugProbe.Touch($"RunBatch start: in={inputDir}, out={outputDir}");
            Directory.CreateDirectory(outputDir);

            var steps = Directory.EnumerateFiles(inputDir, "*.step")
                                 .Concat(Directory.EnumerateFiles(inputDir, "*.stp"))
                                 .ToList();

            _batchQueue = new Queue<string>(steps);
            _batchOutDir = outputDir;

            if (_batchTimer != null)
            {
                _batchTimer.Stop();
                _batchTimer.Dispose();
                _batchTimer = null;
            }

            _batchTimer = new Timer();
            _batchTimer.Interval = 50; // small delay; we do heavy work inside Tick
            _batchTimer.Tick += (s, e) =>
            {
                _batchTimer.Stop(); // stop while processing one part

                try
                {
                    if (_batchQueue.Count == 0)
                    {
                        DebugProbe.Touch("RunBatch done");
                        _batchTimer.Dispose();
                        _batchTimer = null;
                        return;
                    }

                    var step = _batchQueue.Dequeue();
                    DebugProbe.Touch("Processing " + Path.GetFileName(step));

                    ProcessOneStep(step, _batchOutDir); // <— we'll implement next

                    DebugProbe.Touch("Finished " + Path.GetFileName(step));
                }
                catch (Exception ex)
                {
                    DebugProbe.Touch("RunBatch EXCEPTION: " + ex);
                    // Continue to next file
                }
                finally
                {
                    if (_batchTimer != null)
                        _batchTimer.Start();
                }
            };

            _batchTimer.Start();
        }

        private static void ProcessOneStep(string stepPath, string outputDir)
        {
            FeatureCAM.FMDocument doc = null;

            try
            {
                DebugProbe.Touch("Opening STEP via Documents.Open: " + stepPath);
                fc.Documents.Open(stepPath);

                doc = (FeatureCAM.FMDocument)fc.ActiveDocument;

                DebugProbe.Touch("Cast to FMDocument: SUCCESS");

                // 🔹 COMPUTE METRICS HERE
                var metrics = ComputeCncMetrics(doc, outputDir);

                // 🔹 WRITE JSON
                string jsonPath = Path.Combine(
                    outputDir,
                    Path.GetFileNameWithoutExtension(stepPath) + ".json");
                var metric = FCToCAMplete.ComputeCncMetrics(doc, outputDir);
                File.WriteAllText(jsonPath, metric.ToJson());


                DebugProbe.Touch("Wrote " + jsonPath);
            }
            catch (Exception ex)
            {
                DebugProbe.Touch("ProcessOneStep EXCEPTION: " + ex);
            }
            finally
            {
                if (doc != null)
                {
                    try
                    {
                        doc.Close(false);
                        DebugProbe.Touch("Doc.Close(false) worked");
                    }
                    catch { }
                }
            }
        }

        internal static class DebugProbe
        {
            public static void Touch(string tag)
            {
                var p = Path.Combine(
                    Environment.GetFolderPath(Environment.SpecialFolder.Desktop),
                    "FeatureCAM_Addin_Probe.txt");

                File.AppendAllText(p, $"{DateTime.Now:O}  {tag}\r\n");
            }
        }

        public static void OnConnect(object obj, tagFMAddInFlags flags)
        {
            DebugProbe.Touch("OnConnect start");

            try
            {
                fc = GetFeatureCamApp(obj);
                DebugProbe.Touch("fc assigned: " + (fc != null));

                if (fc == null)
                {
                    DebugProbe.Touch("ERROR: Could not acquire FeatureCAM.Application");
                    // Don't call InitializeVariables or CommandBars if fc is null
                    return;
                }

                InitializeVariables();

                // Defer showing UI
                var t = new System.Windows.Forms.Timer();
                t.Interval = 800;
                t.Tick += (s, e) =>
                {
                    t.Stop();
                    t.Dispose();
                    try
                    {
                        if (main_form == null || main_form.IsDisposed)
                            main_form = new UI();
                        main_form.Show();
                        main_form.BringToFront();
                        main_form.Activate();
                        DebugProbe.Touch("UI shown (deferred)");
                    }
                    catch (Exception ex)
                    {
                        DebugProbe.Touch("EXCEPTION showing UI (deferred): " + ex);
                    }
                };
                t.Start();

                // Only now is it safe to call fc.CommandBars...
                // (and only create the button once)
            }
            catch (Exception ex)
            {
                DebugProbe.Touch("EXCEPTION in OnConnect: " + ex);
            }
        }
        
        public static void OnDisConnect(tagFMAddInFlags flags)
        {
            fc.CommandBars.DeleteButton("Utilities", "FeatureCAMToCAMplete");
            if (flags == tagFMAddInFlags.eAIF_DisConnectUserUnLoad)
            {
            }
        }
        private static object GetProp(object comObj, string propName)
        {
            try
            {
                return comObj.GetType().InvokeMember(
                    propName,
                    BindingFlags.GetProperty,
                    null,
                    comObj,
                    null);
            }
            catch
            {
                return null;
            }
        }

        public static void FeatureCAMToCAMplete()
        {
            Variables.is_export_project = true;
            Variables.output_dirpath = "";

            FeatureCAM.FMDocument doc;
            doc = (FeatureCAM.FMDocument)fc.ActiveDocument;
            InitializeVariables();

            if (Variables.doc == null)
            {
                MessageBox.Show("No files are open", Variables.prog_name);
                return;
            }

            // helper function to force a single instance of plugin form
            if (main_form != null)
            {
                main_form.BringToFront();
            }
            else
            {
                main_form = new UI();
                main_form.Show();
                main_form.TopLevel = true;
                main_form.TopMost = true;
                System.Windows.Forms.Application.Run(main_form);
            }
        }

        private static void InitializeVariables()
        {
            FeatureCAM.FMSetup setup = null;

            if (fc != null)
                Variables.docObj = fc.ActiveDocument;
            if (Variables.doc == null)
            {
                Variables.prev_doc_name = "";
                Variables.output_dirpath = "";
            }
            else
            {
                Variables.stock = Variables.doc.Stock;
                Variables.setup_names = new List<string>();
                for (int i = 1; i <= Variables.doc.Setups.Count; i++)
                {
                    setup = (FeatureCAM.FMSetup)Variables.doc.Setups.Item(i);
                    if (setup != null)
                    {
                        Variables.setup_names.Add(setup.Name);
                        /* Have to subtract 1 b/c setups are 1-based and combobox values are 0-based */
                        if (Variables.doc.ActiveSetup.Name == setup.Name)
                            Variables.selected_setup_id = i - 1;
                    }
                }
                Variables.orig_single_stock = Variables.stock.SingleProgramWithProgramStop;

				Variables.clamps = new List<SolidInfo>();
				foreach (FeatureCAM.FMSolid solid in Variables.doc.Solids)
					Variables.clamps.Add(new SolidInfo(solid, solid.UseAsClamp));

					if (Variables.stock.IndexType != FeatureCAM.tagFMIndexType.eIT_None)
					Variables.output_dirpath = Path.Combine(Variables.doc.path, Variables.doc.PartName);
				else
					Variables.output_dirpath = Path.Combine(Variables.doc.path, Variables.doc.PartName) + "_" + Variables.setup_names[Variables.selected_setup_id];

				Variables.doc.ActiveSetup.GetMachineSimLocation(out Variables.offset_x, out Variables.offset_y, out Variables.offset_z);

                Variables.prev_doc_name = Variables.doc.Name;
            }
        }

        public static void Convert() {
            SetupInfo setup_info;
            string tool_info = "";
            bool all_setups_milling;

            try
            {
                Variables.docObj = fc.ActiveDocument;

                if (Variables.stock.IndexType == tagFMIndexType.eIT_None)
                {
                    Variables.doc.Setups.Item(Variables.selected_setup_id + 1).Activate();
                    Variables.stock.SingleProgramWithProgramStop = false;
                }
                /* Verify that file is open */
                if (Variables.doc == null)
                {
                    MessageBox.Show("No files are open", Variables.prog_name);
                    return;
                }

                all_setups_milling = true;
                foreach (FeatureCAM.FMSetup setup in Variables.doc.Setups)
                {
                    if (setup.Type != FeatureCAM.tagFMSetupType.eST_Milling)
                        all_setups_milling = false;
                }
                if (!all_setups_milling)
                {
                    MessageBox.Show("The addin doesn't support non-Milling setups yet. Cannot continue since your part has non-Milling setups.", Variables.prog_name);
                    return;
                }

                Variables.is_single_program =
                    (
                        (Variables.stock.IndexType == FeatureCAM.tagFMIndexType.eIT_None &&
                            Variables.stock.SingleProgramWithProgramStop)
                        ||
                        (Variables.stock.IndexType != FeatureCAM.tagFMIndexType.eIT_None &&
                            !Variables.stock.ToolDominant &&
                            Variables.stock.SingleProgram)
                        ||
                        (Variables.stock.IndexType != FeatureCAM.tagFMIndexType.eIT_None &&
                            Variables.stock.ToolDominant)
                    );


                Directory.CreateDirectory(Variables.output_dirpath);

                /* Initialize setup information (set enabled/disabled) */
                setup_info = null;
                foreach (FeatureCAM.FMSetup setup in Variables.doc.Setups)
                {
                    if (!Variables.is_single_program)
                    {
                        setup_info = new SetupInfo(setup);
                        if (Variables.setups_info == null) Variables.setups_info = new List<SetupInfo>();
                        Variables.setups_info.Add(setup_info);
                        if (setup_info.enabled && setup_info.num_features > 0)
                            Variables.num_enabled_setups++;
                    }
                    else
                    {
                        if (setup_info == null)
                        {
                            setup_info = new SetupInfo(setup);
                            setup_info.name = Variables.doc.Name + "_combined_setup";
                            if (Variables.setups_info == null) Variables.setups_info = new List<SetupInfo>();
                            Variables.setups_info.Add(setup_info);
                        }
                        else
                        {
                            if (Variables.setups_info.Count == 1)
                                Variables.setups_info[0].ucss.Add(new UCS(setup.ucs));
                            else
                                MessageBox.Show("Program failed to generate setup information.\n" + Variables.output_msg, Variables.prog_name);
                        }
                    }
                }

                if (!SaveNCCode())
                {
                    MessageBox.Show("SaveNCCode returned false");
                    return;
                }

                tool_info = ToolsToXmlFile();
                File.WriteAllText(Path.Combine(Variables.output_dirpath, "tools.tdb"), tool_info);

                ExportStock();
                ExportPartSolid();
                ExportClamps();

                CreateCAMpleteProject();
            }
            catch (Exception Ex)
            {
                MessageBox.Show(Ex.Message, "From Convert");
            }
            finally
            {
                Variables.stock.SingleProgramWithProgramStop = Variables.orig_single_stock;
                Variables.Cleanup();
            }
        }
        
        private static void CreateCAMpleteProject()
        {
            string fpath = Path.Combine(Variables.output_dirpath,
                                        Path.GetFileNameWithoutExtension(Variables.doc.Name) + ".proj");

            StringBuilder fcontent = new StringBuilder();

            string units = (Variables.doc.Metric ? "MM" : "INCHES");

            fcontent.AppendLine("<PROJECTCONFIG>");
            fcontent.AppendLine(Lib.tab + "<SOURCE>");
            fcontent.AppendLine(Lib.double_tab + "<CAMSYSTEM>FeatureCAM</CAMSYSTEM>");
            fcontent.AppendFormat(Lib.double_tab + "<VERSION>{0}</VERSION>\n", fc.Version);
            fcontent.AppendLine(Lib.tab + "</SOURCE>");
            fcontent.AppendFormat(Lib.tab + "<NAME>{0}</NAME>\n", Path.GetFileNameWithoutExtension(Variables.doc.Name));
            fcontent.AppendLine(Lib.tab + "<TOOLING>");
            fcontent.AppendFormat(Lib.double_tab + "<TOOLLIBRARY LOADER=\"FEATURECAM_XML_TOOLING\">{0}</TOOLLIBRARY>\n", ".\\tools.tdb");
            fcontent.AppendLine(Lib.tab + "</TOOLING>");
            fcontent.AppendLine(Lib.tab + "<TOOLPATHS>");
            if (!Variables.is_single_program && Variables.stock.IndexType != tagFMIndexType.eIT_None)
            {
                for (int i = 0; i < Variables.setups_info.Count; i++)
                {
                    if (Variables.setups_info[i].enabled && Variables.setups_info[i].num_features > 0)
                    {
                        fcontent.AppendFormat(Lib.double_tab + "<TOOLPATH LOADER=\"FEATURECAM_ACL\" UNITS=\"{0}\">{1}</TOOLPATH>\n", units, Variables.setups_info[i].nc_fpath.Replace(Variables.output_dirpath, "."));
                    }
                }
            }
            else
            {
                fcontent.AppendFormat(Lib.double_tab + "<TOOLPATH LOADER=\"FEATURECAM_ACL\" UNITS=\"{0}\">{1}</TOOLPATH>\n", units, Variables.setups_info[0].nc_fpath.Replace(Variables.output_dirpath, "."));
            }
            fcontent.AppendLine(Lib.tab + "</TOOLPATHS>");
            fcontent.AppendFormat(Lib.tab + "<OFFSETS>\n");
            fcontent.AppendFormat(Lib.double_tab + "<OFFSET TYPE=\"GCODETOPALLETSHIFT\" X=\"{0}\" Y=\"{1}\" Z=\"{2}\" UNITS=\"INCH\"></OFFSET>\n",
                                Variables.offset_x, Variables.offset_y, Variables.offset_z, units);
            fcontent.AppendFormat(Lib.tab + "</OFFSETS>\n");
            fcontent.AppendLine(Lib.tab + "<PARTINFO>");
            fcontent.AppendFormat(Lib.double_tab + "<STOCK LOADER=\"GENERIC_STL\" UNITS=\"{0}\">{1}</STOCK>\n",
                                  units, Variables.stock_fpath.Replace(Variables.output_dirpath, "."));
            if (Variables.clamp_fpaths != null)
                if (Variables.clamp_fpaths.Count > 0)
                {
                    foreach (string clamp_fpath in Variables.clamp_fpaths)
                        fcontent.AppendFormat(Lib.double_tab + "<FIXTURE LOADER=\"GENERIC_STL\" UNITS=\"{0}\">{1}</FIXTURE>\n",
                                              units, clamp_fpath.Replace(Variables.output_dirpath, "."));
                }
            if (Variables.is_export_part)
            {
                if (!String.IsNullOrEmpty(Variables.part_fpath))

                    fcontent.AppendFormat(Lib.double_tab + "<TARGETMODEL LOADER=\"GENERIC_STL\" UNITS=\"{0}\">{1}</TARGETMODEL>\n",
                                      units, Variables.part_fpath.Replace(Variables.output_dirpath, "."));
            }
            fcontent.AppendLine(Lib.tab + "</PARTINFO>");
            fcontent.AppendLine("</PROJECTCONFIG>");

            File.WriteAllText(fpath, fcontent.ToString());
        }

        private static void ExportPartSolid()
        {
            string err_str;

            if (!Variables.is_export_part_now) return;

            FeatureCAM.FMSolid solid = (FeatureCAM.FMSolid)Variables.doc.Solids.Item(Variables.part_solid_name);
            if (solid != null)
            {
                Variables.part_fpath = Path.Combine(Variables.output_dirpath, solid.Name + "_part.stl");
                solid.ExportToSTL(Variables.part_fpath, out err_str);
                if (err_str == "" || err_str == null)
                    Variables.output_msg += Variables.part_fpath + "\n";
                else
                    MessageBox.Show("Error occured while exporting part solid to .stl file: \n'" + err_str + "'", Variables.prog_name);
            }
        }

        private static void ExportStock()
        {
            string err_msg = "";

            try
            {
                Variables.stock_fpath = Path.Combine(Variables.output_dirpath, "stock.stl");

                Variables.stock.ExportToSTL(Variables.stock_fpath, out err_msg);
                if (String.IsNullOrEmpty(err_msg))
                    Variables.output_msg += Variables.stock_fpath + "\n";
                else
                    MessageBox.Show("Error occurred while exporting stock to .stl file: \n'" + err_msg + "'", Variables.prog_name);

            }
            catch (Exception Ex)
            {
                MessageBox.Show("Exception occurred: " + Ex.Message, Variables.prog_name);
            }
        }

        private static void ExportClamps()
        {
            string err_msg = "";
            string fpath;
            try
            {
                Variables.clamp_fpaths = null;
                if (Variables.clamps == null)
                {
                    return;
                }
                if (Variables.clamps.Count == 0)
                {
                    return;
                }
                foreach (SolidInfo clamp in Variables.clamps)
                {
                    if (clamp.is_export)
                    {
                        fpath = Path.Combine(Variables.output_dirpath, clamp.solid.Name + "_clamp.stl");
                        clamp.solid.ExportToSTL(fpath, out err_msg);
                        if (String.IsNullOrEmpty(err_msg))
                            Variables.output_msg += Variables.stock_fpath + "\n";
                        else
                            MessageBox.Show("Error occurred while exporting stock to .stl file: \n'" + err_msg + "'", Variables.prog_name);
                        if (Variables.clamp_fpaths == null) Variables.clamp_fpaths = new List<string>();
                        Variables.clamp_fpaths.Add(fpath);
                    }
                }
            }
            catch (Exception Ex)
            {
                MessageBox.Show("Exception occurred: " + Ex.Message, Variables.prog_name);
            }
        }


        private static bool SaveNCCode()
        {
            int nc_files_num, doc_files_num, macro_files_num;
            object doc_file_names, macro_file_names, nc_file_names;
            string err_msg;
            bool is_op_error;
            int nc_file_id;

            is_op_error = false;
            foreach (FeatureCAM.FMOperation op in Variables.doc.Operations)
                if (!String.IsNullOrEmpty(op.Errors.Trim())) is_op_error = true;

            if (is_op_error)
            {
                MessageBox.Show("Cannot export data to CAMplete: there are errors in the document and nc code cannot be generated.", Variables.prog_name);
                return false;
            }

            /* Set correct output units */
            fc.PostOptionsMill.SetIsInchOutputUnits(!Variables.doc.Metric);
            fc.PostOptionsTurn.SetIsInchOutputUnits(!Variables.doc.Metric);
            fc.PostOptionsWire.SetIsInchOutputUnits(!Variables.doc.Metric);

            /* If part is non-indexed, we can only generate CAMplete report for the active setup */
            if (Variables.stock.IndexType == tagFMIndexType.eIT_None)
                Variables.doc.SaveNC("nc_program.acl", Variables.output_dirpath, false,
                                    FeatureCAM.tagFMSaveNCFileType.eNCFT_NCCode, false, out err_msg,
                                    out nc_files_num, out nc_file_names, out doc_files_num, out doc_file_names,
                                    out macro_files_num, out macro_file_names);
            else if (Variables.stock.SingleProgram) /* We'll have NC code for one file */
                Variables.doc.SaveNC("nc_program.acl", Variables.output_dirpath, false,
                                    FeatureCAM.tagFMSaveNCFileType.eNCFT_NCCode, false, out err_msg,
                                    out nc_files_num, out nc_file_names, out doc_files_num, out doc_file_names,
                                    out macro_files_num, out macro_file_names);
            else
            {
                foreach (FMSetup setup in Variables.doc.Setups)
                {
                    if (setup.Enabled)
                    {
                        setup.Activate();
                        Variables.doc.SimToolpath(false);
                    }
                }
                Variables.doc.SaveNC("nc_program.acl", Variables.output_dirpath, false,
                    FeatureCAM.tagFMSaveNCFileType.eNCFT_NCCode, false, out err_msg,
                    out nc_files_num, out nc_file_names, out doc_files_num, out doc_file_names,
                    out macro_files_num, out macro_file_names);
            }

            if (!Variables.is_single_program && Variables.stock.IndexType != tagFMIndexType.eIT_None)
            {
                if ((int)nc_files_num == Variables.num_enabled_setups)
                {
                    nc_file_id = 1;
                    for (int i = 0; i < Variables.setups_info.Count; i++)
                    {
                        if (Variables.setups_info[i].enabled && Variables.setups_info[i].num_features > 0)
                        {
                            Variables.setups_info[i].nc_fpath = (string)(((Array)nc_file_names).GetValue(nc_file_id));
                            Variables.output_msg += Variables.setups_info[i].nc_fpath + "\n";
                            nc_file_id++;
                        }
                    }
                }
            }
            else
            {
                if ((int)nc_files_num == 1)
                {
                    for (int i = 0; i < Variables.setups_info.Count; i++)
                        Variables.setups_info[i].nc_fpath = (string)(((Array)nc_file_names).GetValue(1));
                    Variables.output_msg += Variables.setups_info[0].nc_fpath + "\n";
                }
            }
            return true;
        }

        private static string ToolsToXmlFile()
        {
            List<string> setup_tools = null;
            List<string> partline_features = null;
            string tool_info = "";
            int setup_num;
            string setup_tool_list;
            FeatureCAM.FMToolMap2 toolmap;

            setup_num = 0;
            setup_tool_list = "";
            partline_features = GetAllFeaturesUsingPartLineComp();
            foreach (FeatureCAM.FMSetup setup in Variables.doc.Setups)
            {
                if (setup_tools == null) setup_tools = new List<string>();
                foreach (FeatureCAM.FMFeature feat in setup.Features)
                    foreach (FeatureCAM.FMOperation op in feat.Operations)
                        setup_tool_list += op.Tool.Name + ";";
                setup_tools.Add(setup_tool_list);
                setup_num++;
            }

            Variables.unsupported_tool_names = "";
            Variables.doc.InvalidateAll();
            /* If we need to create separate tls file for each setup, write tools for each setup to a separate file */
            if (!Variables.is_single_program)
            {
                for (int si = 1; si <= Variables.doc.Setups.Count; si++)
                {
                    tool_info = "";
                    for (int i = 1; i <= Variables.doc.ToolMaps.Count; i++)
                    {
                        toolmap = Variables.doc.ToolMaps.Item(i);
                        if (setup_tools[si - 1].IndexOf(toolmap.Tool.Name + ";") >= 0)
                        {
                            tool_info +=
                                 Tool.ToString(toolmap, partline_features) + Environment.NewLine;
                        }
                    }
                }
            }
            else
            {
                for (int i = 1; i <= Variables.doc.ToolMaps.Count; i++)
                {
                    toolmap = Variables.doc.ToolMaps.Item(i);
                    tool_info +=
                        Tool.ToString(toolmap, partline_features) + Environment.NewLine;
                }
            }

            if (!String.IsNullOrEmpty(Variables.unsupported_tool_names))
                MessageBox.Show("Warning: Tools info was exported, but information for following tool(s) " +
                                Variables.unsupported_tool_names +
                                " was not exported completely, because the tool group(s) are unsupported by this addin.", Variables.prog_name);

            return "<TOOLDB VER=\"2\">" + Lib.EOL +
                        tool_info + //already had end of line at the end
                    "</TOOLDB>";
        }

        private static List<string> GetAllFeaturesUsingPartLineComp()
        {
            List<string> partline_features = null;
            FMLinearPattern lin_ptrn;
            FMRadialPattern rad_ptrn;
            FMPointListPattern ptlist_ptrn;
            FMRectPattern rect_ptrn;
            string pattern_from;
            tagFMFeatureType featureType;
            object attr_val = null;

            foreach (FMFeature feat in Variables.doc.Features)
            {
                pattern_from = "";
                attr_val = null;
                feat.GetFeatureType(out featureType);
                switch (featureType)
                {
                    case tagFMFeatureType.eFT_LinearPattern:
                        lin_ptrn = (FMLinearPattern)feat;
                        attr_val = lin_ptrn.Object.get_Attribute(tagFMAttributeId.eAID_PartLineProgram, null);
                        pattern_from = lin_ptrn.Object.Name;
                        break;
                    case tagFMFeatureType.eFT_RadialPattern:
                        rad_ptrn = (FMRadialPattern)feat;
                        attr_val = rad_ptrn.Object.get_Attribute(tagFMAttributeId.eAID_PartLineProgram, null);
                        pattern_from = rad_ptrn.Object.Name;
                        break;
                    case tagFMFeatureType.eFT_RectPattern:
                        rect_ptrn = (FMRectPattern)feat;
                        attr_val = rect_ptrn.Object.get_Attribute(tagFMAttributeId.eAID_PartLineProgram, null);
                        pattern_from = rect_ptrn.Object.Name;
                        break;
                    case tagFMFeatureType.eFT_PtListPattern:
                        ptlist_ptrn = (FMPointListPattern)feat;
                        attr_val = ptlist_ptrn.Object.get_Attribute(tagFMAttributeId.eAID_PartLineProgram, null);
                        pattern_from = ptlist_ptrn.Object.Name;
                        break;
                    default:
                        attr_val = feat.get_Attribute(tagFMAttributeId.eAID_PartLineProgram, null);
                        break;
                }
                if (attr_val != null)
                {
                    if (System.Convert.ToBoolean(attr_val) == false)
                    {
                        if (partline_features == null) partline_features = new List<string>();
                        partline_features.Add(feat.Name);
                        if (!String.IsNullOrEmpty(pattern_from))
                            partline_features.Add(pattern_from);
                    }
                }
            }
            return partline_features;
        }
    }
}