import argparse
import calendar
import glob
import re

HEADING_UNDERLINE = {
  "version" : "-",
  "section" : "^"
}

SECTIONS_ICONS = {
  "Added"     : "🌱",
  "Changed"   : "⚠️",
  "Deprecated": "🧊",
  "Removed"   : "❌",
  "Fixed"     : "🔧",
}

IGNORE_EMPTY_SECTIONS = True

documented_features = []

def list_documented_features(directory):
  """
  Get list of documented features
  """
  global documented_features

  pattern_ref = r'^\.\. _.*:'
  regex_ref = re.compile(pattern_ref)

  for path_to_file in glob.glob(directory + '/**/*.rst', recursive=True):
    try:
      with open(path_to_file, 'r', encoding='utf-8') as file:
        for line in file:
          match = regex_ref.search(line)
          if match:
            # remove ".. _" prefix and ":" suffix
            documented_features.append(match.group(0).strip()[4:-1])
    except Exception as e:
      # print(f"Error when processing {path_to_file}: {e}")
      pass




def highlight(matchobj):
  """
  Highlight words using custom RST role
  """
  word = matchobj.group(0)
  tail = ""
  if word[-1] == ")" and "(" not in word:
    tail = word[ -1]
    word = word[:-1]
  if word[4:] in documented_features:
    return f":ref:`{word}<{word[4:]}>`"
  else:
    return f":pdmkw:`{word}`{tail}"


def parse_markdown(file_in):
  """
  Convert (pseudo) Markdown file into a structured ChangeLog (dict)
  """

  # Load Markdown changelog as lines of text
  with open(file_in, "r") as f:
    txt_in = f.read().split("\n")

  # Parse lines
  changelog = dict()
  for line in txt_in:
    if line.startswith("## "):
      # Start new version
      str_version = line[line.find("[")+1:line.find("]")]
      if line.find(" - ") < 0:
        date = None
      else:
        date = line[line.find(" - ")+3:].split("-")

      assert(str_version not in changelog)
      changelog[str_version] = dict()
      version = changelog[str_version]

      version["date"]     = date
      version["sections"] = dict()


    elif line.startswith("### "):
      # Start new section
      str_section = line[4:]
      assert(str_section not in version["sections"])

      version["sections"][str_section] = list()
      section = version["sections"][str_section]

      started_list = False


    elif len(line) > 0:
      # New line of section content

      if started_list:
        prev_indent = indent
      else:
        prev_indent = 0

      if line.lstrip().startswith("-"):
        # this is a list item
        started_list = True
        indent = line.find("-")
        if indent != prev_indent:
          section.append("") # add blank line to mark different indentation level

        section.append(line)
      else:
        if started_list:
          # a list is open, append current line to previous list item
          section[-1] += " " + line.lstrip()
        else:
          # plain line (not in a list)
          section.append(line)

  return changelog


def render_rst(changelog, ignore_empty_sections=False, ignore_before_version=None):
  """
  Convert changelog (dict) into ReStructured Text
  """
  pattern = r'PDM_\w+\.*\(*\w+\)*'

  txt_out = []
  indent  = ""
  #  Versions
  for i_version, v in enumerate(changelog):

    if ignore_before_version is not None:
      if v < str(ignore_before_version):
        continue

    date = changelog[v]["date"]
    if date is None: # Version without a date (development version)
      heading = v
    else: # Version with a date
      year, month, day = [int(x) for x in date]
      heading = f"Version {v} ({calendar.month_name[month]} {year})"
    txt_out.append("")
    if i_version == 0:
      txt_out.append(heading)
      txt_out.append(HEADING_UNDERLINE["version"] * len(heading))
    else:
      if i_version == 1:
        txt_out.append("|") # add larger space between versions
        title = "Earlier versions"
        txt_out.append("")
        txt_out.append(title)
        txt_out.append(HEADING_UNDERLINE["version"] * len(title))
        txt_out.append("")
      txt_out.append(f".. dropdown:: {heading}")
      indent = "  "

    #  Sections
    i_section = 0
    for s in changelog[v]["sections"]:

      if ignore_empty_sections and len(changelog[v]["sections"][s]) == 0:
        # Skip empty section
        continue

      if i_version > 0 and i_section > 0:
        txt_out.append("")
        txt_out.append(indent + "|") # add larger space between sections

      txt_out.append("")
      len_under = len(s)
      section = s

      decoration = "**" if i_version > 0 else ""

      # Add appropriate icon
      if section in SECTIONS_ICONS:
        section = decoration + SECTIONS_ICONS[s] + " " + section + decoration
        len_under += 3 + 2*len(decoration)
      txt_out.append(indent + section)
      if i_version == 0:
        txt_out.append(indent + HEADING_UNDERLINE["section"] * len_under)
      txt_out.append("")

      #  Lines
      for l in changelog[v]["sections"][s]:
        # Highlight PDM_* words
        line = re.sub(pattern, highlight, l, flags=re.IGNORECASE)

        txt_out.append(indent + line)

      i_section += 1


  return "\n".join(txt_out)



if __name__ == "__main__":
  # Parse command line args
  parser = argparse.ArgumentParser()

  parser.add_argument("-i",   "--input",                 type=str, default="../ChangeLog")
  parser.add_argument("-p",   "--path",                  type=str, default="../doc/sphinx/source")
  parser.add_argument("-o",   "--output",                type=str, default="../doc/sphinx/source/changelog.rst")
  parser.add_argument("-ibv", "--ignore_before_version", type=str, default=None)
  parser.add_argument("-v", "--verbose", action="store_true")

  args = parser.parse_args()

  # Get list of documented features
  list_documented_features(args.path)
  if args.verbose:
    print(f"{documented_features=}")


  # Load input file
  changelog = parse_markdown(args.input)

  # Print structured ChangeLog
  if args.verbose:
    for v in changelog:
      print(f"\n{v}")
      for s in changelog[v]["sections"]:
        print(f"  {s}")
        for l in changelog[v]["sections"][s]:
          print(f"    {l}")

  # Generate RST str
  rst_txt = render_rst(changelog, IGNORE_EMPTY_SECTIONS, args.ignore_before_version)

  # Print RST output
  if args.verbose:
    print(rst_txt)

  # Write RST output
  with open(args.output, "w") as out:
    out.write(".. _change_log:\n")
    out.write("\n")
    out.write("Release notes\n")
    out.write("#############\n")
    out.write(rst_txt)
