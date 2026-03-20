"""
- download the constellation data from Jonathan McDowell's page, 
- parse it, and 
- convert into a JSON format compatible with our constellations.json structure. 

Output file:  "constellations_mcDowell_<DATE>.json"
"""



mcDowellURL = "https://planet4589.org/space/con/conlist.html"

# download the page
import requests
from pathlib import Path

local_html = Path("conlist.html")
if local_html.exists():
	html_content = local_html.read_text(encoding="utf-8")
else:
	response = requests.get(mcDowellURL, timeout=30)
	response.raise_for_status()
	html_content = response.text
	local_html.write_text(html_content, encoding="utf-8")

# parse the HTML
from bs4 import BeautifulSoup
soup = BeautifulSoup(html_content, 'html.parser')

# for each constellation in the main table, eg 
# "SG1: Starlink Constellation, Modified Gen 1( Apr 2020Filing)"
# "KP2: Kuiper Gen2", etc, extract the data, 
# including the subshells, eg for SG1, the subshells are:
#  "Group 1", "Shell 1-1", "Shell 1-2"...

import json
import re
from datetime import datetime, timezone


def _parse_header(text):
	match = re.match(r"^\s*([A-Za-z0-9]+)\s*:\s*(.*)$", text)
	if not match:
		return None
	code = match.group(1).strip()
	rest = match.group(2).strip()
	name = rest
	comment = ""
	if "(" in rest:
		name = rest.split("(", 1)[0].strip()
		comment = rest[len(name):].strip()
		if comment.startswith("(") and comment.endswith(")"):
			comment = comment[1:-1].strip()
	return code, name, comment


def _strip_suffix(value, suffix):
	if value.endswith(suffix):
		return value[: -len(suffix)].strip()
	return value.strip()


def _parse_row_cells(cells, target_len=15):
	if len(cells) < 5:
		return None
	cols = []
	for i in range(len(cells) - 1):
		cols.append(_strip_suffix(cells[i], cells[i + 1]))
	cols.append(cells[-1].strip())
	if len(cols) < target_len:
		cols.extend([""] * (target_len - len(cols)))
	return cols


def _parse_number(value, as_int=False):
	value = value.strip()
	if not value or value == "-":
		return None
	try:
		num = float(value)
	except ValueError:
		return None
	if as_int:
		return int(round(num))
	return num


def _slugify(value):
	return re.sub(r"[^a-z0-9]+", "", value.lower())


def _make_shell_key(code, layer, subshell, layer_name, existing):
	code_slug = _slugify(code)
	layer_slug = _slugify(layer) or "layer"
	subshell_slug = _slugify(subshell) or "0"
	name_slug = _slugify(layer_name) or "shell"
	base = f"c{code_slug}{layer_slug}{subshell_slug}{name_slug}"
	key = base
	counter = 2
	while key in existing:
		key = f"{base}{counter}"
		counter += 1
	return key


def _shell_label(code, layer_name, layer, subshell):
	if layer_name:
		return f"{code}.{layer_name}"
	return f"{code}.L{layer}S{subshell}"


def _is_total_row(cols):
	if not cols:
		return False
	for value in cols[:5]:
		if not value:
			continue
		upper = value.upper()
		if upper == "TOTAL" or upper.startswith("TOTAL ") or " TOTAL" in upper:
			return True
	return False


tables = soup.find_all('table')
target_table = None
for table in tables:
	headers = [th.get_text(strip=True) for th in table.find_all('th')]
	if any('Constellation' in header for header in headers):
		target_table = table
		break

if target_table is None:
	raise RuntimeError("Could not find constellation table on McDowell page.")

constellations = []
current_constellation = None
shell_keys = set()

for row in target_table.find_all('tr'):
	cells = [cell.get_text(' ', strip=True) for cell in row.find_all(['th', 'td'])]
	if not cells:
		continue

	header = _parse_header(cells[0])
	if header:
		code, name, comment = header
		current_constellation = {
			"name": name,
			"comment": comment,
			"code": code,
			"shells": {},
			"total": None
		}
		constellations.append(current_constellation)
		shell_keys = set()
		continue

	if current_constellation is None:
		continue

	cols = _parse_row_cells(cells)
	if cols is None:
		continue

	constellation_name = cols[0]
	if not constellation_name:
		continue

	layer = cols[2]
	subshell = cols[3]
	layer_name = cols[4]
	if _is_total_row(cols):
		current_constellation["total"] = {
			"layer": layer,
			"subshell": subshell,
			"layer_name": layer_name,
			"alt": _parse_number(cols[5]),
			"inc": _parse_number(cols[6]),
			"nPlane": _parse_number(cols[7], as_int=True),
			"nSat": _parse_number(cols[8], as_int=True),
			"sats_on_station_subshell": _parse_number(cols[9], as_int=True),
			"sats_on_station_layer": _parse_number(cols[10], as_int=True),
			"sats_off_station": _parse_number(cols[11], as_int=True),
			"sats_down": _parse_number(cols[12], as_int=True),
			"total_sats_launched": _parse_number(cols[13], as_int=True),
			"total_sats_planned": _parse_number(cols[14], as_int=True)
		}
		continue

	alt = _parse_number(cols[5])
	inc = _parse_number(cols[6])
	n_plane = _parse_number(cols[7], as_int=True)
	n_sat = _parse_number(cols[8], as_int=True)
	total_planned = _parse_number(cols[14], as_int=True)

	if total_planned is None:
		for idx in (9, 8, 7):
			total_planned = _parse_number(cols[idx], as_int=True)
			if total_planned is not None:
				break

	if total_planned is None:
		continue
	if total_planned == 0:
		continue

	if n_plane is None:
		n_plane = 1
	if n_sat is None:
		n_sat = total_planned

	shell_key = _make_shell_key(current_constellation["code"], layer, subshell, layer_name, shell_keys)
	shell_keys.add(shell_key)
	current_constellation["shells"][shell_key] = {
		"label": _shell_label(current_constellation["code"], layer_name, layer, subshell),
		"totSat": total_planned,
		"nPlane": n_plane,
		"nSat": n_sat,
		"inc": inc,
		"alt": alt
	}





# convert in json format, following the format of constellations.json, eg
# [
#   {
#     "name" : "Starlink Gen-1 (old)",
#     "code" : "SL1old",
#     "comment" : "Starlink-1 former architecture, deprecated",
#     "shells" : {
#         "csl1a1" : {
#             "label" :"SL1 B0",
#             "totSat": 1584,
#             "nPlane": 72,
#             "nSat"  : 32,
#             "inc"   : 53.0,
#             "alt"   : 550.0},
#         "csl1a3" : {
#             "label" :"SL1 B3",
#             "totSat": 75 ,
#             "nPlane": 10,
#             "nSat"  : 38,
#             "inc"   : 81.0,
#             "alt"   : 560.0},
#         "csl1a4" : {
#             "label" :"SL1 B4",
#             "totSat": 450,
#             "nPlane": 12,
#             "nSat"  : 38,
#             "inc"   : 70.0,
#             "alt"   : 570.0},
# etc...


# save the json to a file, eg "constellations_mcDowell_<DATE>.json"

output_name = f"constellations_mcDowell_{datetime.now(timezone.utc).strftime('%Y%m%d')}.json"
with open(output_name, "w", encoding="utf-8") as handle:
	json.dump(constellations, handle, indent=4)

print(f"Saved {len(constellations)} constellations to {output_name}")