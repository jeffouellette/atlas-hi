# ATLAS Latex

ATLAS LaTeX class, style files and templates to typeset notes and papers.
See ChangeLog or Git log for history of changes.

*Responsible:* Ian Brock (Ian.Brock@cern.ch)

Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration

------

## Included files

The following template main files exist:

  - `atlas-paper.tex`: ATLAS paper draft (including CONF and PUB notes)
  - `atlas-note.tex`: ATLAS note
  - `atlas-book.tex`: Long ATLAS document,  such as a TDR
  - `atlas-draft-cover.tex`: Make a standalone cover for an ATLAS draft
  - `atlas-preprint-cover.tex`: Make a standalone cover for an ATLAS CERN preprint
  - `atlas-hepdata-main.tex`:	A front page for material destined for HEPData
  
The ATLAS document class (`atlasdoc.cls`) and style files can be found in 
the latex directory. The following main style files exist:

  - `atlasbiblatex.sty`: Reference style adjustments for biblatex
  - `atlascover.sty`: Make a cover (CONF note, CERN preprint, ATLAS draft)
  - `atlascontribute.sty`: List of contributors (and authors) for a document
  - `atlaspackage.sty`: Standard packages used in ATLAS documents
  - `atlasphysics.sty`: Useful definitions. This file simply inputs others.

Options can be used to specify which should be included.

Documentation can be found in the doc directory.

  - `atlas-latex.pdf`:		Guide to the use of the ATLAS document templates and styles
  - `atlas-bibtex.pdf`:		Guide to references and BibTeX in ATLAS
  - `atlas-physics.pdf`:	Symbols defined in atlasphysics.sty
  - `atlas-tables.pdf`:	Some guidelines and help for tables  

More detailed information about the package can be found under:

<https://twiki.cern.ch/twiki/bin/view/AtlasProtected/PubComLaTeX>

## How to use

The general idea is that, for each document, this package should be cloned into a new directory.
It is assumed that all style files are in a directory latex, which is a subdirectory of 
the one in which the main document sits.

The latex subdirectory can of course be a link to a central style directory.

You can use the `\ATLASLATEXPATH` variable to specify an arbitrary directory.  

To make a new paper/CONF note/PUB note draft give the command:

	make newpaper [BASENAME=mydocument] [TEXLIVE=YYYY]

To make a new ATLAS note give the command:

	make newnote [BASENAME=mydocument] [TEXLIVE=YYYY]

`make new` is an alias for `make newpaper`.

The TeX Live version is set to 2016 by default.
This version number should be fine for newer versions of TeX Live 
and also for an up-to-date MikTeX 2.9 installation. The command `make help` gives you a bit more assistance on which make targets exist.

To add the cover pages for a paper/CONF note/PUB note when circulating it
to the ATLAS collaboration, add the option `coverpage` to the `\documentclass`.

If you want to use the templates for documents that are stored in CERN GitLab,
but are not inside PO-GitLab and hence should not make use of the PO-GitLab CI tools,
you should delete the file: `.gitlab-ci.yml`.
For PO-GitLab documents, the Git repository is in a subdirectory of: https://gitlab.cern.ch/groups/atlas-physics-office/subgroups.


### Running on lxplus
The most common FAQ I get is why atlaslatex does not just compile "out of the box"?
If you are running on `lxplus` for it to work, you MUST set your PATH correctly as follows:

	export PATH=/afs/cern.ch/sw/XML/TL2016/bin/x86_64-linux:$PATH

in order to use TeX Live 2016.
TeX Live 2016 is the current version used in the Physics Office Continuous Integration.

A users guide to the templates can be found in `doc/atlas_latex.pdf`. You can produce
this document yourself (and thus test that your LaTeX setup is working)
by giving the commands:

	cd doc/atlas_latex
	make

or  

	cd doc/atlas_latex
	pdflatex atlas_latex
	biber    atlas_latex
	pdflatex atlas_latex
	pdflatex atlas_latex

Three other make targets are:

  - `make clean`: Cleans up intermediate files
  - `make cleanpdf`: Remove output pdf file
  - `make cleanall`: Also cleans up output pdf file
