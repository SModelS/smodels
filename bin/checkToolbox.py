#!/usr/bin/python

if __name__ == "__main__":
  import set_path, argparse, types
  from Tools import ToolBox
  argparser = argparse.ArgumentParser(description='simple script to check if the tools are installed and compiled')
  ToolBox.toolBox.checkInstallation(colors=True )
