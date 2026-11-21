# Local installation and development

Run these commands in Mac Terminal, starting from your existing repository
directory containing `python-package/`. Python 3.10 or later is required.

```bash
python3 --version
python3 -m venv python-package/.venv
source python-package/.venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e "./python-package[dev]"
python -m pytest python-package/tests
python python-package/examples/quick_start.py
```

The package's `.gitignore` excludes this `.venv/` environment. Keep the
environment on your computer; do not upload it. Do not upload `dist/`, `build/`,
`__pycache__/`, `.pytest_cache/` or `*.egg-info/` directories.

The editable installation makes changes to the package's source available
in this environment without reinstalling. Restart an existing Python or
Jupyter session after changing imported source. For ordinary use, an
installation without `-e` is sufficient.

## Build installable distributions

```bash
cd python-package
python -m build
```

This creates a wheel and source archive under `dist/`. A wheel can be
installed with `python -m pip install path/to/file.whl`. Building does not
publish the package to PyPI or GitHub.

## Update GitHub

From the repository root, first review the changes:

```bash
git status --short
git diff -- python-package
git add python-package
git commit -m "Add Python package for anytime-valid t-tests"
git pull --rebase origin main
git push origin main
```

If a rebase reports conflicts, resolve them before continuing; do not use
a force push to bypass them. If you have other uncommitted work, commit
or otherwise preserve it before the pull.

The existing `r-package/` and `python-package/` folders coexist. Users
install the implementation they need through its subdirectory.
