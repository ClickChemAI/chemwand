# ChemWand

**Tell Claude what to draw. Get an editable ChemDraw file.**

ChemWand is an [MCP](https://modelcontextprotocol.io) server that lets Claude (or any MCP client)
produce publication-style chemistry figures as **CDXML** files — the native, fully editable
format of ChemDraw. You describe the figure in plain language; Claude turns it into SMILES and a
layout description; ChemWand lays it out to ACS 1996 drawing settings, checks that every structure
reads back to the molecule that was asked for (InChI round-trip), and writes the file.

It runs entirely on your own computer. No account, no API key, no server. (The only network call
is the optional PubChem lookup by compound name.)

```
You:     Draw a Suzuki coupling: 4-bromoanisole + phenylboronic acid,
         Pd(PPh3)4, K2CO3, toluene/H2O, 90 °C → 4-methoxybiphenyl, 92%.
Claude:  → draw_linear_scheme(...)  → the .cdxml file opens in ChemDraw
```

## What it can draw

| tool | what it makes |
|---|---|
| `draw_molecule` | a single structure from SMILES |
| `draw_linear_scheme` | reaction schemes, one row or several, with conditions above/below arrows |
| `draw_snaking_scheme` | long multi-step routes folded into equal-width rows |
| `draw_scope_table` | substrate-scope tables (the product grid with yields, ee, dr …) |
| `draw_catalytic_cycle` | catalytic cycles and mechanism rings |
| `draw_radial_figure` | one starting material → many products |
| `merge_figures` | combine figures already drawn into one file |
| `stereo_ingest`, `stereo_lint`, `stereo_propagate`, `stereo_gallery` | read, audit and propagate stereochemistry across a scheme |
| `molecular_properties` | formula, MW, exact mass, and other numbers a paper reports |
| `lookup_structure` | the structure of a named compound, from PubChem |
| `reaction_amounts` | how much of each reagent to weigh out at a given scale |
| `check_disconnection` | atom bookkeeping for one retrosynthetic disconnection |
| `write_procedure` | an ACS-SI-style experimental paragraph for one step |

Every drawing is checked before it is written: a structure whose drawn form does not read back to
the requested molecule is reported, not silently drawn wrong.

## Requirements

- **macOS** or **Windows** (developed and tested on macOS 15; Windows support is newer and less tested)
- **Python 3.10+** (tested on 3.13)
- An MCP client — e.g. [Claude Code](https://claude.com/claude-code) or Claude Desktop
- **ChemDraw** to view and edit the output. Without it the `.cdxml` files are still written;
  any program that reads CDXML can open them.

## Install

### The easy way: let Claude do it

Download this repository (green **Code** button → *Download ZIP*, or `git clone`), open the folder
in Claude Code, and say:

> Install ChemWand from this folder as a local MCP server, following the README.

### By hand

```bash
git clone https://github.com/ClickChemAI/chemwand.git
cd chemwand
python3 -m venv venv
venv/bin/pip install -r requirements.txt          # Windows: venv\Scripts\pip install -r requirements.txt
```

`requirements.txt` pins the exact versions ChemWand was tested with. If a pin has no wheel for
your Python, install the direct dependencies instead: `pip install mcp rdkit numpy pillow uharfbuzz`.
`uharfbuzz` matters — without it, text widths are approximated and layouts shift slightly.

Check it loads (expect `16 tools`):

```bash
venv/bin/python -c "import asyncio, chemwand_mcp as m; print(len(asyncio.run(m.mcp.list_tools())), 'tools')"
```

### Register it with your client

**Claude Code** (use absolute paths):

```bash
# macOS
claude mcp add --scope user chemwand -- /path/to/chemwand/venv/bin/python /path/to/chemwand/chemwand_mcp.py
```

```powershell
# Windows (PYTHONUTF8=1 is required)
claude mcp add --scope user -e PYTHONUTF8=1 chemwand -- C:\path\to\chemwand\venv\Scripts\python.exe C:\path\to\chemwand\chemwand_mcp.py
```

**Claude Desktop** — add to `claude_desktop_config.json`
(macOS: `~/Library/Application Support/Claude/`, Windows: `%APPDATA%\Claude\`), then fully quit and
reopen Claude Desktop. On Windows use `venv\\Scripts\\python.exe` and doubled backslashes.

```json
{
  "mcpServers": {
    "chemwand": {
      "command": "/path/to/chemwand/venv/bin/python",
      "args": ["/path/to/chemwand/chemwand_mcp.py"],
      "env": {"PYTHONUTF8": "1"}
    }
  }
}
```

Quit Claude Desktop before editing this file: it rewrites the file while running and can discard
your change.

### Optional settings (environment variables)

| variable | default | |
|---|---|---|
| `CHEMWAND_FIG_DIR` | `~/chemwand_figures` | where figures are saved |
| `CHEMDRAW_APP` | `ChemDraw 23.1.2.app` | macOS only: your ChemDraw's application name, **including `.app`**. Figures still open through the default `.cdxml` handler if this is wrong; only the optional clean-up step needs it. |
| `CHEMDRAW_FONT` | auto-detected Arial | path to an Arial `.ttf`, if it is not found automatically |

## Use

Restart your client after registering, then just ask:

- *Draw ibuprofen.*
- *Draw this three-step route as one scheme: …*
- *Make a substrate-scope table for this reaction; R = Me, Et, iPr, Ph, 4-MeOC6H4, with these yields: …*
- *Draw the catalytic cycle of a Buchwald–Hartwig amination.*
- *How much of each reagent do I need for this step at 0.5 mmol?*

Figures are saved to `~/chemwand_figures/` and opened automatically. Use a new file name for each
version — reopening a file under the same name can leave ChemDraw showing the old window.

## Notes and limitations

- **Abbreviations.** A hand-written table of 41 common groups (Boc, Cbz, Fmoc, Troc, Alloc, Bn, Ph,
  OAc, OMe, OTBS, OTIPS, OTBDPS, OTMS, OPMB, OMOM, OSEM, OTHP, NHBoc, CF3, NO2, …). Groups outside
  it are drawn out in full.
- The layout engine is tuned to ACS 1996 settings (14.4 pt bonds, 10 pt Arial).
- Windows: file opening, file locking and font lookup have Windows code paths, but most testing was
  on macOS. Bug reports welcome.

## License

MIT — see [LICENSE](LICENSE). Third-party notices are in [THIRD_PARTY_NOTICES.txt](THIRD_PARTY_NOTICES.txt).

ChemDraw is a trademark of Revvity, Inc. ChemWand is an independent project and is not affiliated
with, endorsed by, or sponsored by Revvity. It writes CDXML files; it does not include or modify any
ChemDraw software.
