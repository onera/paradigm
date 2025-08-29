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
      date        = line[line.find(" - ")+3:].split("-")

      assert(str_version not in changelog.keys())
      changelog[str_version] = dict()
      version = changelog[str_version]

      version["date"]     = date
      version["sections"] = dict()


    elif line.startswith("### "):
      # Start new section
      str_section = line[4:]
      assert(str_section not in version["sections"].keys())

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


def render_rst(changelog, ignore_empty_sections=False):
  """
  Convert changelog (dict) into ReStructured Text
  """
  pattern = r'PDM_\w+\.*\(*\w+\)*'

  txt_out = []
  #  Versions
  for v in changelog.keys():
    txt_out.append("")
    if len(txt_out) > 1: txt_out.append("|") # add larger space between versions
    txt_out.append("")
    year, month, day = [int(x) for x in changelog[v]["date"]]
    heading = f"Version {v} ({calendar.month_name[month]} {year})"
    txt_out.append(heading)
    txt_out.append(HEADING_UNDERLINE["version"] * len(heading))

    #  Sections
    for s in changelog[v]["sections"].keys():

      if ignore_empty_sections and len(changelog[v]["sections"][s]) == 0:
        # Skip empty section
        continue

      txt_out.append("")
      len_under = len(s)
      section = s
      # Add appropriate icon
      if section in SECTIONS_ICONS:
        section = SECTIONS_ICONS[s] + " " + section
        len_under += 3
      txt_out.append(section)
      txt_out.append(HEADING_UNDERLINE["section"] * len_under)

      #  Lines
      for l in changelog[v]["sections"][s]:
        # Highlight PDM_* words
        line = re.sub(pattern, highlight, l, flags=re.IGNORECASE)

        txt_out.append(line)

  return "\n".join(txt_out)



if __name__ == "__main__":
  # Parse command line args
  parser = argparse.ArgumentParser()

  parser.add_argument("-i", "--input",   type=str, default="../ChangeLog")
  parser.add_argument("-p", "--path",    type=str, default="../doc/sphinx/source")
  parser.add_argument("-o", "--output",  type=str, default="../doc/sphinx/source/changelog.rst")
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
    for v in changelog.keys():
      print(f"\n{v}")
      for s in changelog[v]["sections"].keys():
        print(f"  {s}")
        for l in changelog[v]["sections"][s]:
          print(f"    {l}")

  # Generate RST str
  rst_txt = render_rst(changelog, IGNORE_EMPTY_SECTIONS)

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
