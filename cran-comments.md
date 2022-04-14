Dear CRAN

this is a resubmission of the R package MLFS. I have received the first comments
from CRAN with the instructions and suggestions for improvement. I would like to
thank you for your time and effort in reviewing my package. I have implemented 
all the suggestions and am now submitting the updated and corrected version. 
Please read my responses and described actions below. I have also performed all
the necessary tests listed at the end of this message.

Best,
Jernej


CRAN cooments

* 1. If there are references describing the methods in your package, please 
add these in the description field of your DESCRIPTION file in the form
authors (year) <doi:...>
authors (year) <arXiv:...>
authors (year, ISBN:...)
or if those are not available: <https:...>
with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for 
auto-linking.
(If you want to add a title as well please put it in quotes: "Title")

  * Author: Two references are now added in the DESCRIPTION file. 


* 2. Please add \value to .Rd files regarding exported methods and explain 
the functions results in the documentation. Please write about the 
structure of the output (class) and also what the output means. (If a 
function does not return a value, please document that too, e.g. 
\value{No return value, called for side effects} or similar).
Missing Rd-tags in up to 25 .Rd files, e.g.:
      add_stand_variables_halfPeriod.Rd: \value
      add_stand_variables.Rd: \value
      BAI_prediction_halfPeriod.Rd: \value
      BAI_prediction.Rd: \value
      calculate_BAL_halfPeriod.Rd: \value
      calculate_BAL.Rd: \value
      ...

  * Author: I have now described the outputs of all functions. I used @return tag. 

* 3. \ dontrun{} should only be used if the example really cannot be executed 
(e.g. because of missing additional software, missing API keys, ...) by 
the user. That's why wrapping examples in \ dontrun{} adds the comment 
("# Not run:") as a warning for the user.
Does not seem necessary.
Please unwrap the examples if they are executable in < 5 sec, or replace 
\dontrun{} with \donttest{}.

  * Author: The example which was previously wrapped in \ dontrun{}, is now wrapped in \ donttest{}. The reason for wrapping is because it is not executable in < 5 sec.



* 4. If you use a package which is only needed in examples, please list it in 
'Suggests' and wrap these examples in if(requireNamespace("pkgname")){} 
instead.

  * Author: In the originally submitted version, I had an R package that was only needed in one example, lmfor. However, the R package lmfor was recently removed from CRAN. Therefore, I disabled the functionality related to the R package lmfor and also removed the example. There are no other R package dependencies that would be used only in the examples.

* 5. You write information messages to the console that cannot be easily 
suppressed.
It is more R like to generate objects that can be used to extract the 
information a user is interested in, and then print() that object.
Instead of print()/cat() rather use message()/warning()  or 
if(verbose)cat(..) (or maybe stop()) if you really have to write text to 
the console.
(except for print, summary, interactive functions)
Please fix and resubmit.
  * Author: I use now message() instead of print(). 

##  Resubmission submission
* This is a re-submission of the package MLFS after the first CRAN comments.

## Test environments
* local OS X install, R 4.1.1

* rhub Windows Server 2022 (https://builder.r-hub.io/status/original/MLFS_0.3.9.tar.gz-528bf27f02674d8396b073f911a356ad)
* rhub Ubuntu (https://builder.r-hub.io/status/original/MLFS_0.3.9.tar.gz-c4868e3e6b79456fb27d50facaf55c53)
* rhub Fedora Linux (https://builder.r-hub.io/status/original/MLFS_0.3.9.tar.gz-bf76c3f96b08429bb163f25b0e661ae7)

* win-check oldrelease (https://win-builder.r-project.org/Ot2b4WJk8daC/00check.log)
* win-check release (https://win-builder.r-project.org/kcU9eqG7W7h1/00check.log)
* win-check devel (https://win-builder.r-project.org/SsmZKi7IB5I6/00check.log)

## R CMD check results
There were 0 ERRORs, 0 WARNINGs and 0 NOTEs

## Downstream dependencies
* Will be checked in all future versions
