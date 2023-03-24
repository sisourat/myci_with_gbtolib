import sys
import os
import readline
import argparse
from gooey import Gooey, GooeyParser

ormas_pyt_path='/home/nico/Workspace/GeneralFanoCICode/PythonScript/ormas_like.py'
running = True

@Gooey(optional_cols=2, program_name="GUI Python script to general CI input files")
def main():

  settings_msg = 'The script will help you to set up the calculations'
  parser = GooeyParser(description=settings_msg)
  parser.add_argument('--verbose', help='be verbose', dest='verbose', action='store_true', default=False)


# setup CI calculations
  subs = parser.add_subparsers(help='commands', dest='command')
  ci_parser = subs.add_parser('CI_calculations', help='')

  ci_parser.add_argument('ormas', metavar='Name of the ormas python script?  ', help='', widget='FileChooser', default=ormas_pyt_path)
  ci_parser.add_argument('fxml', metavar='Name of the XML file?  ', help='', widget='FileChooser')
  ci_parser.add_argument('inputname', metavar='Name of the CI file(s) that will be generated?  ', help='', default='input')
  ci_parser.add_argument('fshell', metavar='Name of the script that will be generated?  ', help='', widget='FileChooser')


  ci_parser.add_argument('--nci', metavar='For multiple CI  ', help='Number of CI?', default=1)
  ci_parser.add_argument('--imin', metavar='', help='First index?', default=-1)

# setup Total decay widths calculations
  tgam_parser = subs.add_parser('Total_decay_widths', help='')
  tgam_parser.add_argument('anxml', metavar='How many XML files do you have?  ', help='Enter an integer')

  args = parser.parse_args()

# CI calculations

  scriptname=args.fshell
  fxml=args.fxml
  nci=int(args.nci)
  imin=int(args.imin)
  inputname=args.inputname
  ormaspath=args.ormas

  print ""
  print "A shell script will be generated"
  fshell = open(scriptname,"w")
  ninput = -1

  if(nci==1):
    print >> fshell, "python ", ormaspath, fxml, "> a"
    print >> fshell, "cat header.txt a > ", inputname
  else:
    for iorb in range(imin,imin+nci):
      ninput+=1
      print >> fshell, 'sed "s/X/'+str(iorb)+'/g"', fxml, "> a"
      print >> fshell, "python ", ormaspath," a > b"
      print >> fshell, "cat header.txt b > ", inputname+str(ninput)

  print >> fshell, "rm a b header.txt"
  cmd="chmod u+x "+scriptname
  os.system(cmd)
  print "Dont forget to run ./"+scriptname
  print "Run myCI < input_myci"

  print args.fxml
  print args.fshell
  print args.nci
  print args.imin


if __name__ == '__main__':
        main()

