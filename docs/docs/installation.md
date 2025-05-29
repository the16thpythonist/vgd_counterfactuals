# Installation

This guide will help you install VGD Counterfactuals and its dependencies.

## Requirements

VGD Counterfactuals requires Python 3.9 or later (up to Python 3.12). Make sure you have a compatible Python version installed:

```bash
python --version
```

## Installation Methods

### Method 1: Install from PyPI (Stable Release)

You can install the stable release of VGD Counterfactuals directly from PyPI:

```bash
pip install vgd_counterfactuals
```

### Method 2: Install from Source (Recommended)

**1. Clone the repository:**
```bash
git clone https://github.com/the16thpythonist/vgd_counterfactuals.git
cd vgd_counterfactuals
```

**2. Create a virtual environment (recommended):**
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

**3. Install the package:**
```bash
pip install -e .
```

## Verify Installation

After installation, verify that VGD Counterfactuals is properly installed:

```bash
python -m vgd_counterfactuals.cli --version
python -m vgd_counterfactuals.cli --help
```

You should see the version information and help text.

## Next Steps

Once you have VGD Counterfactuals installed, head over to the [Get Started](get-started.md) guide to learn how to use the library with a simple example.
