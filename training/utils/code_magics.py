import sys
import glob
import os
import json
from subprocess import run, CalledProcessError
# from tempfile import NamedTemporaryFile as tmpf
from IPython.core import magic_arguments
from IPython.core.magic import cell_magic, line_magic, Magics, magics_class

# from nbformat import read, NO_CONVERT


headers = {
  "c"       : "#include <stdlib.h>\n#include <stdio.h>\n#include <mpi.h>\n\n",
  "fortran" : "",
  "python"  : ""
}


language_extension = {
  "c"      : "c",
  "fortran": "F90",
  "python" : "py"
}

language_compiler = {
  "c"      : "mpicc",
  "fortran": "mpif90",
  "python" : None
}

# temporary files (separate code blocks, not merged files and executables)
# will be placed here:
tmp_dir = "./tmp"



def print_out_err(proc):
  """ Print output and error of a process"""
  stdout = proc.stdout
  stderr = proc.stderr
  if stderr:
    sys.stderr.write(stderr)
    sys.stderr.flush()
  if stdout:
    sys.stdout.write(stdout)
    sys.stdout.flush()



@magics_class
class CodeMagics(Magics):
  """
  Code block with integer identifier (to merge with other code blocks)
  """
  @cell_magic
  @magic_arguments.magic_arguments()
  @magic_arguments.argument('--id', '-i',
                            help='Block identifier (order)',
                            type=int,
                            default=0)
  @magic_arguments.argument('--prefix', '-p',
                            help='File name prefix',
                            type=str,
                            default="tmp")
  @magic_arguments.argument('--language', '-l',
                            help='Language',
                            type=str,
                            choices=["c", "python", "fortran"],
                            default="python")

  def code_block(self, line="", cell=None):
    args = magic_arguments.parse_argstring(self.code_block, line)

    # extension = language_extension[args.language]
    extension = "tmp"
    str_id = "%3.3d" % args.id
    self.filename = f"{tmp_dir}/{args.prefix}_{str_id}.{extension}"

    if not os.path.exists(tmp_dir):
      os.makedirs(tmp_dir)

    # open source file and write cell content
    print(f"code block written in {self.filename}")
    with open(self.filename, "w+b") as out:
      out.write(bytes(cell, "utf-8"))
      out.flush()


  """
  Merge code blocks with given prefix into a new file
  The merged code can also be compiled and run in an MPI environment
  """
  @line_magic
  @magic_arguments.magic_arguments()
  @magic_arguments.argument('--prefix', '-p',
                            help='File name prefix',
                            type=str,
                            default=None)
  @magic_arguments.argument('--language', '-l',
                            help='Language',
                            type=str,
                            choices=["c", "python", "fortran"],
                            default="python")
  @magic_arguments.argument('--n_rank', '-n',
                            help='Number of MPI ranks',
                            type=int,
                            default=0)
  @magic_arguments.argument('--clear', '-c',
                            help='Clear generated files',
                            action="store_true")
  @magic_arguments.argument('--verbose', '-v',
                            help='Print merged source code',
                            action="store_true")
  def merge_code_blocks(self, line):
    # get current environment
    self.env = os.environ.copy() # ?

    args = magic_arguments.parse_argstring(self.merge_code_blocks, line)

    extension = language_extension[args.language]
    source_name = f"{tmp_dir}/{args.prefix}.{extension}"

    # get ordered list of code blocks
    tmp_files = glob.glob(f"{tmp_dir}/{args.prefix}_*.tmp")

    if len(tmp_files) == 0:
      print("No tmp files were found, make sure you run all code_block cells before this one!")
      return

    # merge into a single file
    with open(source_name, "w") as outfile:
      for fname in sorted(tmp_files):
        with open(fname) as infile:
          outfile.write(infile.read())

    if args.verbose:
      with open(source_name, "r") as f:
        print(f"merged :\n{f.read()}")

    run_code = args.n_rank > 0

    # Compile & run
    if run_code:
      # Compile
      if args.language == "python":
        exec_name = source_name
      else:
        exec_name = os.path.join('.', os.path.basename("{:s}.{:s}".format(source_name, "exe")))
        command = []
        command.extend([language_compiler[args.language]])
        command.extend(["-g"])
        command.extend(['-o', os.path.basename(exec_name), source_name])
        # TODO: rajouter ce qui manque pour compiler
        sys.stdout.write(" ".join(command)+"\n")

        try:
          proc = run(command,
                     capture_output=True,
                     shell=False,
                     check=True,
                     env=self.env,
                     encoding='utf8')
          print_out_err(proc)
        except CalledProcessError as e:
          sys.stderr.write(" ".join(e.cmd))
          sys.stderr.write(e.stderr)
          sys.stderr.write(e.stdout)
          return

      # Run
      command = []
      command.extend(["mpirun"])
      command.extend(["-np"])
      command.extend(["%d" % args.n_rank])
      if args.language == "python":
        command.extend(["python3"])
      command.extend([exec_name])

      sys.stdout.write(" ".join(command)+"\n")

      if os.path.isfile(exec_name):
        try:
          proc = run(command,
                     capture_output=True,
                     shell=False,
                     check=True,
                     env=self.env,
                     encoding='utf8')
          print_out_err(proc)
        except CalledProcessError as e:
          sys.stderr.write(" ".join(e.cmd))
          sys.stderr.write(e.stderr)
          sys.stderr.write(e.stdout)
          return
      else:
        sys.stderr.write(f"No executable was produced. It should be {exec_name}")

    # Clear generated files
    if args.clear:
      rm_files = tmp_files + [source_name]
      if exec_name != source_name:
        rm_files += exec_name
      for fname in rm_files:
        if os.path.exists(fname):
          os.remove(fname)

def load_ipython_extension(ipython):
  ipython.register_magics(CodeMagics)
