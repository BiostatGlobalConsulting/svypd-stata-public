# svypd-stata-public
Stata program to estimate survey proportions

This program estimates survey proportions in Stata.  It does not return results to the log window, but provides some stored scalars and matrices for programmers to access.  To view the estimation results, type the command `return list` after running svypd.

The program addresses some current shortcomings of Stata's `svy: proportion` command.  Those shortcomings and how they are addressed are described in the .do file named "Demonstrate advantages of user-written svypd program.do".

Note that svypd.ado was written as part of the World Health Organization's Vaccination Coverage Quality Indicators (VCQI) suite of programs.
More information about VCQI is available at: www.biostatglobal.com/VCQI_resources.html

Contact Dale Rhoda with questions:  Dale.Rhoda@biostatglobal.com or Dale.Rhoda@gmail.com

www.biostatglobal.com

