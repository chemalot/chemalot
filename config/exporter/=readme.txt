Exporter documentation
config file location:
Aestel/config/exporter/*

The export is driven by tree major perl scripts 
  tabExport.pl runs a single sql statment from a configuration file using optional parameters
               and generates a tab separated file with all columns from the statement
  sdfExport.pl  runs a single sql statment from a configuration file using optional parameters
                and generates a sdf file with all columns from the statement
  exporter.pl runs all commands defined in the toolConfig/commands xml element
              and zips all the resulting output files into a singel zip file
              and then sends that zip file to one or more e-mails

config file format:
- main elements:
  1.) <sql> contains sql query to be executed
      use ? for parameters which will be replaced by command line parameters 
      passed to sdfExporter or tabExporter
      
  2.) <toolConfig>
      configuration parameters for exporter.pl
      
      <toolConfig><commands>
      list of csh commands to be executed. some varaibles starting with "$" can 
      be used here cf. help text printed by executing exporter.pl
      
      <email>albertgo,joe@joe.com</email>
      
