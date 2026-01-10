# Onshape ↔ FeatureCAM Import / Export Tooling

# LLM To Do
1. Change parameters in existing features
2. Add new fillets

This repository contains two coordinated components used in an Onshape → FeatureCAM workflow:

1. **Python script (`Export_test.py`)**  
   Authenticates with the **Onshape API** and interfaces with a specific Onshape document identified by a URL.

2. **.NET Framework Plugin (FeatureCAM 2026)**  
   A **FeatureCAM 2026 plugin** that enables import/export functionality inside FeatureCAM.

---

### Python / Onshape
- Python 3.10+ recommended
- pip
- Valid Onshape API credentials
- Access to the target Onshape document

### FeatureCAM Plugin
- FeatureCAM **2026**
- Visual Studio compatible with the .NET Framework used by the plugin
- FeatureCAM SDK / interop assemblies available on the system

---

## Onshape Python Script (`Export_Test.py`)

### Description
`Export_Test.py` connects to the Onshape API using API keys provided via environment variables and operates on a specific Onshape document referenced by URL (e.g., querying data or exporting geometry for downstream use).

### Required Environment Variables
Secrets must **not** be committed to source control.  
Provide them via local ".env" file.

Required:
- `ONSHAPE_ACCESS_KEY`
- `ONSHAPE_SECRET_KEY`

A template is provided in `.env.example`.

---

### Setup

From the repo root:

```bash
cd python
python -m venv .venv

# PowerShell
.\.venv\Scripts\Activate.ps1
# or cmd.exe
# .\.venv\Scripts\activate.bat

pip install -r requirements.txt

python Export_Test.py

