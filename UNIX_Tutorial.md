###Lab Meeting 9/16/2015###

___
####Introduction to UNIX and Command line Bioinformatics####

###There are many ways to set up your Mac, here is how I would do it. (type or copy/paste commands into terminal)###

##Need to install Xcode (Developer Tools that Apple doesn’t install natively but supports):##

`xcode-select --install`

##Install HomeBrew (copy and paste this into terminal):##

```ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)”```

Then setup homebrew:
```
brew doctor
brew tap homebrew/science
```	
##Install some tools (GNU tools, python, CPAN, BLAST, HMMer3):##

	`brew install coreutils python cpanm blast hmmer`

##Install BioPython##
```
pip install biopython
```

Install BioPerl (CPAN is a helper script to install Perl Modules)

`sudo cpanm BioPerl`

Install a good Text Editor, my favorite is Text Wrangler which also has command line tools that are useful
```
https://s3.amazonaws.com/BBSW-download/TextWrangler_4.5.12.dmg
http://pine.barebones.com/files/tw-cmdline-tools-4512.zip
```
