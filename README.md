The open-source datasets and code resources provided in this project are described as follows:

•	Hollow_lumen_structure and Partitioned_lumen_structure contain the raw printing-path data generated 
  after path planning for vascular models with hollow and partitioned lumen cross-sections, respectively.

•	Scaffold contains the raw printing-path data generated from the differentiated scaffold model.

For the scripts:

•	Scaffold_with_MN-CAPP generates multi-nozzle collaborative printing G-code for the differentiated scaffold based on the MN-CAPP strategy.

•	Scaffold_without_MN-CAPP is the control script that generates G-code using a conventional multi-nozzle collaborative printing method.

•	Vascular_printing is a general G-code generation script for multi-nozzle collaborative printing of both vascular models under the MN-CAPP strategy.

To generate the printing G-code for a vascular structure, simply load either the Hollow_lumen_structure or Partitioned_lumen_structure dataset and run the Vascular_printing script.
