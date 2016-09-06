Running Tests
===============

runTest.py can be used to run automated tests for the command line programs.
Run as:

- runTests.py
  to run all tests

- runTests.py dirNames...
  To run only specific tests

Writing Tests
===============

Tests are configured test.xml files contained in subdirectories of the current.
Each test.xml file can contain multiple &lt;test> entries.
A &lt;test> entry can be configured with these attributes:

- in: required
  input file will be supplied on stdin of the command

- out: required
  output file output of command will be written to this file
  stderr will be written to a file with the same name but the file type will
  be modified to .err: t.txt -> t.err

- wdir: optional
  if not given the tests will be preformed in the same directory as the xml file
  all path given including those in other atributes must be relative to 
  wdir or absolute.


A &lt;test> entry must contain the commands to be executed.
The commands are executed with the tcsh shell with the working directory
being wdir.

The test direcotory (containing the test.xml file) can be referenced in 
the command or any of its attributes by using the "$dir" string.

Before execution the ./out directory is removed and recreated in the same 
directory as the test.xml file. This directory is intended to hold
output files such that the
test directory is not cluttered. This directory can be referenced in the 
commands and attributes by using the $outDir string.


A &lt;test> can have &lt;init> and &lt;postprocess> elements containg commands
to be executed before or after the test command.

A &lt;test> entry can include one or &lt;diff> of &lt;diffDir> elments.

A &lt;diff> element will diff a reference file to the output file. e.g.
>   &lt;diff ref="10.ref.tab"/>


A &lt;diffDir> element will diff a refernce directory to another directory:
>    &lt;diffDir refDir='../ref' outDir='.' opts='--exclude="output.*"'/>

