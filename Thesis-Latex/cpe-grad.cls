%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     Filename = cpe-grad.cls v1.0
%     Author   = Peerapon Siripongwutikorn (peerapon@cpe.kmutt.ac.th)
% 
% This is a LaTeX class file for project report and thesis at CPE,
% KMUTT.
% Date: Dec 25, 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\usepackage[thai,english]{babel}
%\usepackage{thswitch}
\usepackage{times}
\usepackage[left=4cm,top=3cm,right=2cm,bottom=2cm]{geometry} 
\usepackage[pdftex]{graphicx} %%% Use the option pdftex instead of
                             %%% dvips if you use pdflatex 
\usepackage{array}
%\usepackage[ps2pdf]{thumbpdf}
% 1. Generage pdf file
% 2. thumbpdf --modes=dvips pdf-file for ps2pdf
%    thumbpdf pdf-file for pdfLaTeX
% 3. regenerate pdf file with -t option
\usepackage{amsfonts,amsmath} %% For math fonts and operations

\renewcommand{\thairmdefault}{angsana}
\renewcommand{\englishrmdefault}{ptm}
%\renewcommand{\textthai}[1]{
%\usefont{LTH}{angsana}{m}{n}\selectfont\small #1}
\renewcommand{\textthai}[1]{\thaitext{\small #1}}

\def\disstitleone{TITLE OF DISSERTATION}      % 1st line of TITLE of dissertation
\def\disstitletwo{}                           % 2nd line of TITLE if needed
\def\disstitlethree{}                         % 3rd line of TITLE if needed
\def\disstitlefour{}                          % 4th line of TITLE if needed
\def\dissauthor{Your Name}                    % Author of dissertation
\def\dissdiplomaone{B.A., Carnegie Mellon University, 1991} % first degree
\def\dissdiplomatwo{M.S., University of Pittsburgh, 1993}   % second degree
\def\dissdiplomathree{}                       % still in school?
\def\dissdiplomafour{}                        % collecting them?
\def\dissdegree{Doctor of Philosophy}         % Name of the degree
\def\dissdegreeabrev{Ph.D.}                   % Abbreviation of the degree
\def\dissyear{1995}                           % Year of submission
\def\dissdefensedate{June 1 2000}             % Date of defense
\def\disscommitteechair{Prof. Francis Crick}  % Name of the Chairman
\def\dissadvisor{Prof. Linus Pauling}    % 2nd committee member
\def\disscommitteetwo{Prof. Julius Axelrod}   % 3rd committee member
\def\disscommitteethree{}   % 4th committee member
\def\disscommitteefour{}                      % 5th committee member

\def\thaiworktype{}
\def\thaidisstitleone{}
\def\thaidisstitletwo{}
\def\thaidisauthor{}
\def\thaidissadvisor{}
\def\thaidissdegree{}
\def\thaidissyear{}
% 
% Include the followng line (uncommented) in your dissertation file
% if your work has more than one appendix. This is what will be written on the 
% cover page preceding the appendices.  By default, it is set to {Appendix}.
% 
% \def\appendixnames{Appendices}              % if more than one appendix
% 
% Your dissertation should be structured the following way:
% \maketitlepage, \makecopyrightpage, \makesignaturepage,
% \abstract{}, \preface{}, \tableofcontents, \listoftables, \listoffigures,
% \chapter{},....., \appendix{},.....,
% \thebibliography (or \bibliography{} with BiBTeX), \theindex
% 
% ************************************************************************

\input report.cls  % most of the stuff is defined here
\typeout{Document Style `CPE thesis/report'.}
\newcommand{\dotrule}[1]{% 
      \parbox[t]{#1}{\dotfill}} 

% **************************************************************************
% Redefinition of some values
% **************************************************************************
\def\ssp{\def\baselinestretch{1.0}\large\normalsize}
\def\dsp{\def\baselinestretch{1.37}\large\normalsize}
%\def\dsp{\def\baselinestretch{1.1}\large\normalsize}
\newif\ifdraftmode \draftmodefalse
\ifdim \overfullrule > 0pt \draftmodetrue \ssp \else \dsp \fi  % draft mode on?
\def\truebaselineskip{\baselineskip}
\def\baselineskipxii{24pt}
\def\baselineskipxi{18.6pt}  % for title and signature page
\def\emptyline#1{\vspace*{#1\baselineskip}}
\newif\iffirstappendix \firstappendixtrue
\newif\iffirstchapter \firstchaptertrue
\def\tocprefix{}
\newif\ifprelude \preludetrue
\newif\ifromanmode \romanmodetrue \pagenumbering{roman}
\def\setarabic{\ifromanmode \romanmodefalse \pagenumbering{arabic} \fi}
\setcounter{secnumdepth}{5}
\topmargin -0.8cm \oddsidemargin 0.5 in
%\textheight 9.0truein \advance\textheight by \topskip \textwidth 6.0truein
\textheight 24.69cm \advance\textheight by \topskip \textwidth 6.0truein
\footskip 30pt \footnotesep 7.7pt
\headheight 12 truept \headsep 25 truept
\onecolumn \parindent 0em
\if@twoside\else\raggedbottom\fi
\widowpenalty=10000 \clubpenalty=10000 \looseness=1

% ***********************************************************************
% Modifications for single spacing and indentation, from Ethan V. Munson's 
% style file for dissertations at the Univ of California
% ***********************************************************************

\def\smallssp{\def\baselinestretch{1.0}\large\small}
\def\scriptsizessp{\def\baselinestretch{1.0}\large\scriptsize}
%\def\tabular{\par\ssp\let\@halignto\@empty\@tabular}
%\expandafter \def\csname tabular*\endcsname #1{\par\ssp\def\@halignto{to#1}\@tabular}
\def\smalltabular{\par\smallssp\let\@halignto\@empty\@tabular}
\def\endsmalltabular{\endtabular}
\expandafter \def\csname smalltabular*\endcsname #1{\par\smallssp\def\@halignto{to#1}\@tabular}
\expandafter \let \csname endsmalltabular*\endcsname = \endtabular
\def\scriptsizetabular{\par\scriptsizessp\let\@halignto\@empty\@tabular}
\def\endscriptsizetabular{\endtabular}
\expandafter \def\csname scriptsizetabular*\endcsname #1{\par\scriptsizessp\def\@halignto{to#1}\@tabular}
\expandafter \let \csname endscriptsizetabular*\endcsname = \endtabular
\long\def\@footnotetext#1{\insert\footins{\ssp\reset@font\footnotesize
    \interlinepenalty\interfootnotelinepenalty\splittopskip\footnotesep
    \splitmaxdepth \dp\strutbox \floatingpenalty \@MM \hsize\columnwidth \@parboxrestore
   \edef\@currentlabel{\csname p@footnote\endcsname\@thefnmark}\@makefntext
    {\rule{\z@}{\footnotesep}\ignorespaces #1\strut}}}
\def\@startsection#1#2#3#4#5#6{\if@noskipsec \leavevmode \fi
   \par \@tempskipa #4\relax \@afterindenttrue
   \ifdim \@tempskipa <\z@ \@tempskipa -\@tempskipa \relax\fi
   \if@nobreak \everypar{}\else
     \addpenalty{\@secpenalty}\addvspace{\@tempskipa}\fi \@ifstar
     {\@ssect{#3}{#4}{#5}{#6}}{\@dblarg{\@sect{#1}{#2}{#3}{#4}{#5}{#6}}}}

% ***********************************************************************
% Redefinition of main dissertation divisions, to comply with KMUTT 
% regulations.
% ***********************************************************************

%%%%%%%%%%%%%%%%%%%%%%%%% pagestyles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\def\ps@plain{\let\@mkboth\markboth\def\@oddfoot{}\def\@evenfoot{}
   \def\@oddhead{\hbox{}\hfil\rm\thepage} \def\@evenhead{\rm\thepage\hfil\hbox{}}}
\def\ps@plainbottom{\let\@mkboth\markboth\def\@oddhead{}\def\@evenhead{}
   \def\@oddfoot{\hfil\rm\thepage\hfil} \def\@evenfoot{\hfil\rm\thepage\hfil}}
\ps@plainbottom

%%%%%%%%%%%%%%%%%%%%%%%%% Chapter font %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%% Chapter heading in the text for \chapter{}
\def\@makechapterhead#1{%\vspace*{.8truein}%
  {\centering{\large\bf \MakeUppercase{\@chapapp{}} \thechapter}\hskip 5\p@  }%\vskip 20\p@
  {\centering\MakeUppercase{\large\bf #1}\par\nobreak\emptyline 2}}
%  {\centering\MakeUppercase{\large\bf #1}\par\nobreak\vskip 10\p@ }}
%\def\@makeschapterhead#1{\vspace*{.8truein}%
%  {\centering{\large\bf #1}\par\nobreak\vskip 40\p@ }}


\def\chapter{\clearpage
%   \ifprelude \pagestyle{plainbottom} \else 
%   \pagestyle{plain}\thispagestyle{plainbottom} \fi
   \global\@topnum\z@ \@afterindenttrue
   \secdef\@chapter\@schapter}
\def\@chapter[#1]#2{\preludefalse\setarabic
%   \pagestyle{plain}\thispagestyle{plainbottom} % added in v2.4
   \pagestyle{plain}\thispagestyle{empty} % No page number in chapter.
   \refstepcounter{chapter}%
   \typeout{\@chapapp\space\thechapter.}% 
%% Put the work CHAPTER in ToC before the first chapter
   \iffirstchapter
       \addcontentsline{toc}{chapterxx}{\mdseries\bf CHAPTER}
       \firstchapterfalse
   \fi
%% Add the chapter heading in ToC
   \addcontentsline{toc}{chapter}{\normalsize\bf\MakeUppercase{\tocprefix {\protect \numberline{\thechapter.}\hskip -0.5em #1}}}
   \chaptermark{#1}%
   %\addtocontents{lof} {\protect\addvspace{0\p@}}
   %\addtocontents{lot} {\protect\addvspace{0\p@}}
   \addtocontents{lof} {\protect\addvspace{0em}}
   \addtocontents{lot} {\protect\addvspace{0em}}
   \pagestyle{plain}\thispagestyle{empty}
   % Must be added after addcontentsline 
   % Otherwise, the page number may not show up..
   \@makechapterhead{#2} \@afterheading}
%% Remove top space from \chapter*{}
%% Chapter heading in the text for \chapter*{}
\def\@makeschapterhead#1{%
  %\vspace*{50\p@}%
  {\parindent \z@ \centering%\raggedright
    \normalfont
    \interlinepenalty\@M
    \Large \bfseries  #1\par\nobreak
    \vskip 20\p@
  }}
\def\@schapter#1{\@makeschapterhead{#1} \@afterheading}

%%%%%%%%%%%%%%% Section heading fonts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\renewcommand\section{\@startsection {section}{1}{\z@}%
                                   {-3.5ex \@plus -1ex \@minus -.2ex}%
				   {0.1ex\@plus.1ex} 
                                   %{2.3ex \@plus.2ex}%
                                   {\normalfont\large\bfseries}}
\renewcommand\subsection{\@startsection{subsection}{2}{\z@}%
                                     {-3.25ex\@plus -1ex \@minus -.2ex}%
                                     {0.1ex \@plus .2ex}%
                                     {\normalfont\large\bfseries}}
\renewcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}%
                                     {-3.25ex\@plus -1ex \@minus -.2ex}%
                                     {0.1ex \@plus .2ex}%
                                     {\normalfont\normalsize\bfseries}}

%%%%%%%%%%%%%%% Table/figure caption font %%%%%%%%%%%%%%%%%%%%%%%%
\setlength\belowcaptionskip{10\p@}
\long\def\@maketablecaption#1#2{
    \global \@minipagefalse
    \hbox to \hsize{\parbox[t]{\hsize}{{\bf #1} #2 \vskip 0.5em}}}
\let\old@floatboxreset\@floatboxreset
\def\table{%
\let\@makecaption\@maketablecaption
%\def\@floatboxreset{%
% \old@floatboxreset
% \@tablesize
%}%
  \let\footnoterule\relax
  \let\contcaption\conttablecaption \@float{table}}
\def\endtable{\end@float}


\long\def\@makecaption#1#2{%
  \vskip\abovecaptionskip
  \sbox\@tempboxa{{\bf #1} #2}%
  \ifdim \wd\@tempboxa >\hsize
    {\bf #1} #2\par
  \else
    \global \@minipagefalse
    \hb@xt@\hsize{\hfil\box\@tempboxa\hfil}%
  \fi
  \vskip\belowcaptionskip}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define how the content lines show up in ToC when using
%%% addcontentsline{}{}{}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Normal chapter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
\renewcommand*\l@chapter[2]{%
  \ifnum \c@tocdepth >\m@ne
    \addpenalty{-\@highpenalty}%
    \vskip 1.0em \@plus\p@
    %\vskip 0.5em \@plus\p@    %% Space between entry in toc
    \setlength\@tempdima{1.5em}%
    \begingroup
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      \leavevmode %\bfseries
      \advance\leftskip\@tempdima
      \hskip -\leftskip
      \bfseries #1\nobreak\hfil \nobreak\hb@xt@\@pnumwidth{\hss #2}\par
      \penalty\@highpenalty
    \endgroup
  \fi}

%% Special lines such as abstract, acknowledgment
%% have no extra spacing between lines.
\newcommand*\l@chapterx[2]{%
  \ifnum \c@tocdepth >\m@ne
    \addpenalty{-\@highpenalty}%
    \vskip 0.2em \@plus\p@    %% No spacing between lines
    \setlength\@tempdima{1.5em}%
    \begingroup
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      \leavevmode %\bfseries
      \advance\leftskip\@tempdima
      \hskip -\leftskip
      \bfseries #1\nobreak\hfil \nobreak\hb@xt@\@pnumwidth{\hss #2}\par
      \penalty\@highpenalty
    \endgroup 
  \fi}

%%% Special toc entry with no number
\newcommand*\l@chapterxx[2]{%
  \ifnum \c@tocdepth >\m@ne
    \addpenalty{-\@highpenalty}%
    \vskip 1em \@plus\p@ % Space between entry
    \setlength\@tempdima{1.5em}%
    \begingroup
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      \leavevmode %\bfseries
      \advance\leftskip\@tempdima
      \hskip -\leftskip
      \bfseries #1\nobreak\hfil \nobreak\hb@xt@\@pnumwidth{\hss}\par
      \penalty\@highpenalty
    \endgroup \vskip -0.8em
  \fi}

%%% Section entries %%%%%%%%%%%%%%
\def\l@section#1#2{\addpenalty{\@secpenalty}  % good place for page break
   %\addvspace{1.0em plus 1pt}  % space above toc entry
   \addvspace{0.2em}  % space above toc entry
   \@tempdima 2em            % width of box holding section number
   \begingroup
     \parindent \z@ \rightskip \@pnumwidth 
     \parfillskip -\@pnumwidth   
     %\sc                  % Small caps
     \leavevmode          % TeX command to enter horizontal mode.       
      \advance\leftskip\@tempdima  %% added 5 Feb 88 to conform to 
      \hskip -\leftskip            %% 25 Jan 88 change to \numberline
      \hskip 2em #1\nobreak\hfil \nobreak\hb@xt@\@pnumwidth{\hss #2}\par
     %\hskip 2em #1\nobreak\hfil \nobreak\hbox to\@pnumwidth{\hss #2}\par
   \endgroup}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\let\l@subsection\l@section
\let\l@subsubsection\l@section

%%% Figure/table entries in list of figures and tables shows up
\def\l@figure#1#2{\addpenalty{\@secpenalty}  % good place for page break
   %\addvspace{1.0em plus 1pt}  % space above toc entry
  \addvspace{0.2em}  % space above toc entry
   \@tempdima 2em            % width of box holding section number
   \begingroup
     \parindent \z@ \rightskip \@pnumwidth 
     \parfillskip -\@pnumwidth   
     %\sc                  % Small caps
     \leavevmode          % TeX command to enter horizontal mode.       
      \advance\leftskip\@tempdima  %% added 5 Feb 88 to conform to 
      \hskip -\leftskip            %% 25 Jan 88 change to \numberline
      \hskip 0em #1\nobreak\hfil \nobreak\hb@xt@\@pnumwidth{\hss #2}\par
     %\hskip 2em #1\nobreak\hfil \nobreak\hbox to\@pnumwidth{\hss #2}\par
   \endgroup}
\let\l@table\l@figure

%%% Appendix content lines
\def\l@apdx#1#2{\addpenalty{\@secpenalty}  % good place for page break
   \addvspace{0.2em}  % space above toc entry
   \@tempdima 2em            % width of box holding section number
   \begingroup
     \parindent \z@ \rightskip \@pnumwidth 
     \parfillskip -\@pnumwidth   
     %\sc                 % Small caps
     \leavevmode          % TeX command to enter horizontal mode.       
      \advance\leftskip\@tempdima  %% added 5 Feb 88 to conform to 
      \hskip -\leftskip            %% 25 Jan 88 change to \numberline
      \hskip 2em #1\nobreak\hfil \nobreak\hb@xt@\@pnumwidth{\hss #2}\par
   \endgroup}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Title page %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\def\schoolname{}
\def\worktype{}  % Type of work..
\def\titleempty{}
\def\maketitlepage{\thispagestyle{empty}\c@page\@ne\baselineskip \baselineskipxii
   \emptyline 1 {\centering 
   \makebox{\MakeUppercase{\disstitleone}} \par
   \ifx \titleempty\disstitletwo \else  
   \makebox{\MakeUppercase{\disstitletwo}} \par
   \fi
   \emptyline 7 \MakeUppercase{\dissauthor} \par \emptyline 7
   %\makebox{\dissdiplomaone} \par \makebox{\dissdiplomatwo} \par
   %\makebox{\dissdiplomathree} \par \makebox{\dissdiplomafour} \par
   A \MakeUppercase{\worktype}\ SUBMITTED IN PARTIAL FULFILLMENT \par
   OF THE REQUIREMENTS FOR \par 
   THE DEGREE OF \MakeUppercase{\dissdegree}\ (COMPUTER ENGINEERING) \par
   \makebox{\MakeUppercase{\schoolname}} \par
%   \makebox{\schoolname} in partial fulfillment \par
%   of the requirements for the degree of \par \dissdegree \par
   KING MONGKUT'S UNIVERSITY OF TECHNOLOGY THONBURI \par 
   \dissyear \par}
   \newpage\baselineskip \truebaselineskip}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Copyright page %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\def\makecopyrightpage{\newpage
   {\centering{\vbox{}\vfill Copyright by \dissauthor\par\dissyear\par \vfill\vbox{}}}\newpage}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Signature page %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\def\thesis{Thesis}
\def\makesignaturepage{
   \newpage\thispagestyle{empty}\baselineskip \baselineskipxi {\centering
   \makebox {\disstitleone} \par 
   \ifx \titleempty\disstitletwo \else  
   \makebox{\disstitletwo} \par
   \fi
   %\makebox {\disstitletwo} \par
   %\makebox {\disstitlethree} \par \makebox {\disstitlefour} \par
   \vfill \dissauthor\ \ \dissdiplomaone \ (Computer Engineering) \par 
   \vfill A \worktype\ Submitted in Partial Fulfillment \par
   of the Requirements for \par
   the Degree of \dissdegree \ (Computer Engineering) \par
   \makebox{\schoolname} \par
   King Mongkut's University of Technology Thonburi \par 
   \dissyear \par}
 \emptyline 1
\begin{flushleft}
\worktype\ Committee \par\par
\begin{tabular*}{\textwidth}{@{}c@{\extracolsep{1.2cm}}l@{}}
   & \\
%%%%% Chairman %%%%%%%%%%%%%%%%%%%%%%%%%
  \dotrule{3.6in} & Chairman \\
  (\disscommitteechair) & \\ 
   & \\
%%%%% Thesis or project advisor %%%%%%%%
\ifx \thesis\worktype
  \dotrule{3.6in} & Member and  \\
  (\dissadvisor) & \worktype\ Advisor\\
   & \\ 
\else
  \dotrule{3.6in} & Member and Special \\
  (\dissadvisor) & Project Study Advisor\\
   & \\ 
\fi
%%%%%%% 3rd committee %%%%%%%%%%%%%%%%%%
  \dotrule{3.6in} & Member \\
  (\disscommitteetwo) & \\
   & \\
%%%%%%% 4th committee %%%%%%%%%%%%%%%%%%
\ifx \titleempty\disscommitteethree \else
  \dotrule{3.6in} & Member \\
  (\disscommitteethree) & \\
  & \\
\fi
%%%%%%%% 5th committee %%%%%%%%%%%%%%%%
\ifx \titleempty\disscommitteefour \else
  \dotrule{3.6in} & Member \\
  (\disscommitteefour) & \\
   & \\
\fi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\end{tabular*}
\end{flushleft}
\begin{center}	
\vspace*{\fill}Copyright reserved
\end{center}	
\newpage\baselineskip \truebaselineskip}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% English Abstract %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\def\abstract{
\newpage
\setcounter{page}{2}
\pagestyle{plain}%\thispagestyle{empty}
\addcontentsline{toc}{chapterx}{\mdseries ENGLISH ABSTRACT}
\begin{flushleft}
\begin{tabular}{@{}lp{0.7\textwidth}}
 \worktype\ Title & \disstitleone\ \\ 
\ifx \titleempty\disstitletwo \else
              & \disstitletwo \\
\fi
 \worktype\ Credits & \disscredit \\
 Candidate & \dissauthor \\
 \worktype\ Advisor & \dissadvisor \\
 Program & \dissdegree \\
 Field of Study & Computer Engineering \\
 Department & Computer Engineering \\
 Faculty & Engineering \\
 B.E. & \thaidissyear 
\end{tabular}
\end{flushleft}
%}
\emptyline 1
\begin{center} 
  Abstract
\end{center}
\par
}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Thai Abstract %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\def\thaiabstract{
\newpage
\addcontentsline{toc}{chapterx}{\mdseries THAI ABSTRACT}
\begin{flushleft}
\begin{tabular}{@{}lp{0.7\textwidth}}
 \textthai{��Ǣ��}\thaiworktype & \thaidisstitleone\ \\ 
\ifx \titleempty\thaidisstitletwo \else
              & \thaidisstitletwo \\
\fi
 \textthai{˹��¡Ե} & \disscredit \\
 \textthai{�����¹} & \thaidissauthor \\
 \textthai{�Ҩ�������֡��} & \thaidissadvisor \\
 \textthai{��ѡ�ٵ�} & \thaidissdegree \\
 \textthai{�Ң��Ԫ�} & \textthai{���ǡ�������������} \\
 \textthai{�Ҥ�Ԫ�} & \textthai{���ǡ�������������} \\
 \textthai{���} & \textthai{���ǡ�����ʵ��} \\
 \textthai{�.�.} & \thaidissyear 
\end{tabular}
\end{flushleft}
%}
\emptyline 1
\begin{center} 
  \textthai{���Ѵ���}
\end{center}
\par
}
\def\endabstract{\newpage}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Acknowledgement %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\def\preface{\chapter*{\large ACKNOWLEDGMENTS}%
\addcontentsline{toc}{chapterx}{\mdseries ACKNOWLEDGMENTS}
\thispagestyle{plain}
}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Content page %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table of contents
\setcounter{tocdepth}{5}  % must be same as secnumdepth
\def\contentsname{CONTENTS}
\def\tableofcontents{\chapter*{\large\MakeUppercase{\contentsname}\@mkboth{}{}%
 \\ \emptyline 1 \hfill\textbf{PAGE}\\[-0.4em]}{\ssp\@starttoc{toc}}
\addcontentsline{toc}{chapterx}{\mdseries CONTENTS}
\pagestyle{plain}\thispagestyle{plain}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% List of table page %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% List of tables
\def\listtablename{LIST OF TABLES}
\def\listoftables{\chapter*{\large \MakeUppercase{\listtablename}\@mkboth{TABLE}{PAGE}\thispagestyle{plain}
\\ \emptyline 1 \normalsize\textbf{TABLE}\hfill\textbf{PAGE}\\[-0.4em] }
   \addcontentsline{toc}{chapterx}{\mdseries LIST OF TABLES} {\ssp\@starttoc{lot}}}
%\def\l@table{}
%\def\l@table{\vskip 1.0em plus\p@\@dottedtocline{1}{1.5em}{2.3em}}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% List of figure page %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\def\listfigurename{LIST OF FIGURES}
\def\listoffigures{\chapter*{\large\MakeUppercase{\listfigurename}\@mkboth{FIGURE}{PAGE}\thispagestyle{plain}\\
   \emptyline 1
   \normalsize\textbf{FIGURE}\hfill\textbf{PAGE}\\[-0.4em]  }
   \addcontentsline{toc}{chapterx}{\mdseries LIST OF FIGURES} {\ssp\@starttoc{lof}}}
\let\l@figure\l@table

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% List of symbol page %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\def\listofsymbols{\chapter*{\large LIST OF SYMBOLS}
%\\ \emptyline 1 \textbf{SYMBOL}\hfill\textbf{UNIT}\\[-1cm]}
\addcontentsline{toc}{chapterx}{\mdseries LIST OF SYMBOLS}\thispagestyle{plain}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% List of vocab page %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\def\listofvocab{\chapter*{\large LIST OF TECHNICAL VOCABULARY AND ABBREVATIONS}%
\addcontentsline{toc}{chapterx}{\mdseries LIST OF TECHNICAL VOCABULARY AND ABBREVATIONS}\thispagestyle{plain}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Appendix pages %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\firstappendixtrue
\def\appendixnames{Appendix} % assume only one appendix
\def\appendixcover#1{%
\newpage\pagestyle{plain}\thispagestyle{empty}
\iffirstappendix %
  \setcounter{chapter}{0} % Reset chapter numbering
  \renewcommand{\thechapter}{\Alph{chapter}}
  \def\@chapapp{\MakeUppercase{\appendixname}}%
  \addcontentsline{toc}{chapter}{\mdseries\bf \MakeUppercase{\appendixnames}}%
  \firstappendixfalse% 
\fi 
   \setcounter{equation}{0} % Reset equation numbering
   \addtocounter{chapter}{1}  % Start the 1st appendix 
   {\centering{\vbox{}\vfill \MakeUppercase{\bf
    APPENDIX \thechapter} \par \MakeUppercase{#1} \vfill\vbox{}}}\newpage}
%   {\centering{\vbox{}\vfill \MakeUppercase{{\normalsize\bf APPENDIX
%\thechapter}} \par \MakeUppercase{#1}} \vfill\vbox{}}}\newpage}

\def\appendix#1{% 
\appendixcover{#1}
\addcontentsline{toc}{apdx}{\thechapter\hskip 1em #1}%
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Bibliography page %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\def\refname{REFERENCES}
\def\referencepage{\chapter*{\large\bf \MakeUppercase{\refname}}%
\addcontentsline{toc}{chapter}{\mdseries\bf \MakeUppercase{\refname}}\thispagestyle{empty}
}
%\def\bibliographycover{\newpage \typeout{\bibname .} \thispagestyle{empty}
%   {\centering{\vbox{}\vfill {\LARGE\bf \bibname} \vfill\vbox{}}}\newpage}
%%% Change the brackets to dot around bib number.%%%%%%%%%%%%%%
\newdimen\bibindent
\bibindent=2em
\def\@biblabel#1{#1.}
\def\thebibliography#1{%\bibliographycover
  %\chapter*{\bibname\@mkboth{}{}}
  %\addcontentsline{toc}{chapter}{\bibname} \ssp
  \referencepage \vspace*{-30\p@}\ssp\pagestyle{plain}	
  \list{\@biblabel{\arabic{enumiv}}}{\settowidth\labelwidth{\@biblabel{#1}}%
    \leftmargin\labelwidth \advance\leftmargin\labelsep
    \usecounter{enumiv} \let\p@enumiv\@empty
    \def\theenumiv{\arabic{enumiv}}}%
    \def\newblock{\hskip .11em plus.33em minus.07em}%
    \sloppy\clubpenalty4000\widowpenalty4000 \sfcode`\.=\@m }
\def\endthebibliography{\ifdraftmode \dsp \fi
  \def\@noitemerr{\@warning{Empty `thebibliography' environment}} \endlist}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Vitae %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\def\vitae{
\clearpage
\chapter*{\large\bf CURRICULUM VITAE}%
\vspace*{-30\p@}\ssp
\addcontentsline{toc}{chapter}{\normalsize\bf CURRICULUM VITAE}
\pagestyle{plain}\thispagestyle{plain}
}

% ************************************************************
% personal stuff

%\input latexsym.sty
\usepackage{latexsym}

\newcounter{definition}[chapter]
\newenvironment{definition}[1]{\refstepcounter{definition}
   \par\hbox{}\par\noindent
   {\bfseries Definition \thechapter.\thedefinition\ #1}
   \par\noindent\ignorespaces}
   {\par\hbox{}\par}

\newenvironment{proof}{
   \par\noindent{\bfseries Proof: }\ignorespaces}
   { $\Box$\par\hbox{}\par}

% *******************************************************************
\endinput

