# ChemWand

**Tell Claude what to draw. Get an editable ChemDraw file.**

[中文安装说明 →](README_zh.md)

ChemWand is an extension for Claude that draws publication-style chemistry figures — structures,
reaction schemes, multi-step routes, substrate-scope tables, catalytic cycles — and saves them as
**CDXML**, the native editable format of ChemDraw. You describe the figure in plain language;
ChemWand lays it out, checks that every structure reads back to the molecule you asked for
(InChI round-trip), saves the file and opens it in ChemDraw.

It runs on your own computer. No account, no API key, no server. (The only network call is the
optional PubChem lookup by compound name.)

## Install (2 minutes)

You need **Claude Desktop** ([download](https://claude.com/download)) and **ChemDraw**. Nothing else —
Claude installs Python and everything ChemWand needs by itself.

1. Download **[chemwand.mcpb](https://github.com/ClickChemAI/chemwand/releases/latest/download/chemwand.mcpb)**.
2. Double-click it. Claude Desktop opens an install page.
3. Click **Install**. The first install downloads Python and the chemistry libraries, so give it a
   minute or two.

That's it. Open a chat and ask:

> Use ChemWand to draw ibuprofen.

**The Install button is grey?** Your Claude account belongs to an organization (company or university)
whose admin only allows approved extensions. Ask your admin to approve ChemWand, or use
[Option B](#option-b-claude-code) below.

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

Things to try:

- *Draw a Suzuki coupling: 4-bromoanisole + phenylboronic acid, Pd(PPh3)4, K2CO3, toluene/H2O,
  90 °C → 4-methoxybiphenyl, 92%.*
- *Make a substrate-scope table for this reaction; R = Me, Et, iPr, Ph, with yields 85, 80, 72, 90%.*
- *Draw the catalytic cycle of a Buchwald–Hartwig amination.*
- *How much of each reagent do I need for this step at 0.5 mmol?*

Figures are saved to `chemwand_figures` in your home folder and open in ChemDraw automatically. Ask
for a new file name for each version — reopening a file under the same name can leave ChemDraw
showing the old window.

## Other ways to install

### Option B: Claude Code

For accounts that cannot install extensions. Download this repository (green **Code** button →
*Download ZIP*), unzip it, open the folder in the **Code** tab of Claude Desktop (or in Claude Code),
and say:

> Install ChemWand from this folder as a local MCP server, following the README. If Python 3.11 or
> newer is missing, install it in my user folder with uv.

Then fully quit and reopen Claude, and draw from the Code tab.

### Option C: by hand

Requires Python 3.11–3.14 (tested on 3.13).

```bash
git clone https://github.com/ClickChemAI/chemwand.git
cd chemwand
python3 -m venv venv
venv/bin/pip install -r requirements.txt          # Windows: venv\Scripts\pip install -r requirements.txt
venv/bin/python -c "import asyncio, chemwand_mcp as m; print(len(asyncio.run(m.mcp.list_tools())), 'tools')"   # expect 16
```

Register it with Claude Code (absolute paths):

```bash
# macOS
claude mcp add --scope user chemwand -- /path/to/chemwand/venv/bin/python /path/to/chemwand/chemwand_mcp.py
```

```powershell
# Windows (PYTHONUTF8=1 is required)
claude mcp add --scope user -e PYTHONUTF8=1 chemwand -- C:\path\to\chemwand\venv\Scripts\python.exe C:\path\to\chemwand\chemwand_mcp.py
```

Optional environment variables: `CHEMWAND_FIG_DIR` (where figures are saved, default
`~/chemwand_figures`), `CHEMDRAW_APP` (macOS: your ChemDraw's application name including `.app`,
used only by the optional clean-up step), `CHEMDRAW_FONT` (path to an Arial `.ttf` if it is not found
automatically).

## Notes and limitations

- **Abbreviations.** A hand-written table of 41 common groups (Boc, Cbz, Fmoc, Troc, Alloc, Bn, Ph,
  OAc, OMe, OTBS, OTIPS, OTBDPS, OTMS, OPMB, OMOM, OSEM, OTHP, NHBoc, CF3, NO2, …). Groups outside
  it are drawn out in full.
- The layout follows ACS 1996 settings (14.4 pt bonds, 10 pt Arial).
- Developed and tested on macOS. Windows has its own code paths for opening files, file locking and
  fonts, but has had little testing — bug reports welcome in
  [Issues](https://github.com/ClickChemAI/chemwand/issues).

## License

MIT — see [LICENSE](LICENSE). Third-party notices are in [THIRD_PARTY_NOTICES.txt](THIRD_PARTY_NOTICES.txt).

ChemDraw is a trademark of Revvity, Inc. ChemWand is an independent project and is not affiliated
with, endorsed by, or sponsored by Revvity. It writes CDXML files; it does not include or modify any
ChemDraw software.
