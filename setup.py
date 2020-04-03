#!/usr/bin/env python3

"""
.. module:: setup
   :synopsis: Setup script for SModelS.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import os
import sys
from setuptools import setup
from setuptools.command.install import install
sys.path.insert ( 0, "./" )
from smodels.installation import version, authors, requirements, resolve_dependencies, fixpermissions
import subprocess

class OverrideInstall(install):

    def run(self):
        #uid, gid = 0, 0
        install.run(self) # calling install.run(self) insures that everything 
                # that happened previously still happens, 
        install_as_user = True
        try:
            import getpass
            if getpass.getuser() == "root":
                install_as_user = False
        except ImportError:
            pass
        enableStupidMacFix=False
        if enableStupidMacFix:
            if "Apple" in sys.version:
                # a wild attempt at fixing a problem with Mac OS X. Somehow
                # setup.py doesnt resolve the requirements!
                try:
                    self.do_egg_install()
                except NameError:
                    pass
        # so the installation does not break! 
        # here we start with doing our overriding and private magic ..
        mode = 0o777
        for filepath in self.get_outputs():
            # if self.install_scripts in filepath:
            if "smodels/lib/" in filepath:
                #print ("Overriding setuptools mode of scripts ...")
                #log.info("Changing ownership of %s to uid:%s gid %s" %
                #         (filepath, uid, gid))
                #os.chown(filepath, uid, gid)
                # print ("Changing permissions of %s to %s" %
                #         ( os.path.dirname ( filepath ), oct(mode)))
                os.chmod( os.path.dirname ( filepath ), mode )
        if install_as_user: ## FIXME doesnt work for system installs
            resolve_dependencies( as_user=install_as_user )
        if not install_as_user:
            fixpermissions()

def read(fname):
    """
    Simple method to read a file (fname) located in the current folder.

    """
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


def listDirectory (dirname):
    """ list all files in directory, skip subdirs """
    if dirname[-1] == "/":
        dirname = dirname[:-1]
    files = os.listdir (dirname)
    ret = []
    for filename in files:
        fullname = dirname + "/" + filename
        extension = os.path.splitext ( filename )[1]
        if os.path.isdir ( fullname ) or \
                extension in [ ".out", ".tgz", ".1" ] or \
                filename in [ "Makefile", "README" ]:
            continue
        ret.append ( fullname )
    return ret

def dataFiles ():
    """
    List all config files and binaries

    """
    ret = []
    ret.append ( ("smodels/", [ "smodels/version", "smodels/COPYING", "smodels/README.rst", "smodels/INSTALLATION.rst" ]) )
    for directory in [ "smodels/share/", "smodels/share/models/",
          "smodels/etc/", "smodels/lib/nllfast/nllfast-1.2/", 
          "smodels/lib/nllfast/nllfast-2.1/", "smodels/lib/nllfast/nllfast-3.1/", 
          "smodels/lib/pythia6/", "smodels/lib/pythia8/" ]:
        ret.append ((directory, listDirectory (directory)))
    for directory in ["inputFiles/slha/", "inputFiles/lhe/" ]:
        ret.append (( "smodels/"+directory, listDirectory (directory)))

    return ret


def compile():
    """
    Compile external tools by calling make

    """

    if len(sys.argv) < 2:
        return
    needs_build = False
    for i in sys.argv[1:]:
        if i in [ "build", "build_ext", "build_clib", "install", 
                  "install_lib", "bdist", "bdist_rpm", "bdist_dumb", 
                  "bdist_wininst", "bdist_wheel", "develop"]:
            needs_build = True
    if not needs_build:
        return
    subprocess.call(["make", "-C", "smodels/lib" ])

# compile() ## not needed anymore as we perform compilation-on-demand now

setup(
    name = "smodels",
    version = version(),
    author = authors(),
    author_email="smodels-developers@lists.oeaw.ac.at ",
    entry_points = {
            'console_scripts': ['smodels-config=smodels.installation:main',
                           'runSModelS.py=smodels.tools.runSModelS:main',
                           'smodelsTools.py=smodels.tools.smodelsTools:main' ]
    },
    install_requires=requirements(),
    data_files=dataFiles() ,
    description=("A tool for interpreting simplified-model results from the "
                   "LHC"),
    cmdclass={'install': OverrideInstall},
    license="GPLv3",
    # use_2to3 = True,
    keywords=("simplified models LHC BSM theories interpretation "
                "supersymmetry UEDs"),
    url="https://smodels.github.io/",
    packages=['smodels',
              'smodels.theory',
              'smodels.tools',
              'smodels.experiment'],
    include_package_data = True,
    test_suite='test',
    long_description=read('smodels/README.rst'),
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    ]
)
