
/*
---------------------------------------------------------------
ARWEN v1.2.3 Dean Laslett
---------------------------------------------------------------

    ARWEN
    Detects metazoan mitochondrial tRNA genes in nucleotide sequences
    Copyright (C) 2007-2008 Dean Laslett

    Please, report bugs and suggestions of improvements to the authors

    E-mail: Björn Canbäck: bcanback@acgt.se
            Dean Laslett: gaiaquark@gmail.com

    Version 1.2.3  May 19th, 2008.


    Please reference the following paper if you use this
    program as part of any published research.

    Laslett, D. and Canback, B. (2007)
    ARWEN: a program to detect tRNA genes in
    metazoan mitochondrial nucleotide sequences.
    Bioinformatics, 24(2); 172-175.

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; version 2 of the License, (see below).

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.


               GNU GENERAL PUBLIC LICENSE
               Version 2, June 1991

 Copyright (C) 1989, 1991 Free Software Foundation, Inc.
 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 Everyone is permitted to copy and distribute verbatim copies
 of this license document, but changing it is not allowed.

       Preamble

  The licenses for most software are designed to take away your
freedom to share and change it.  By contrast, the GNU General Public
License is intended to guarantee your freedom to share and change free
software--to make sure the software is free for all its users.  This
General Public License applies to most of the Free Software
Foundation's software and to any other program whose authors commit to
using it.  (Some other Free Software Foundation software is covered by
the GNU Library General Public License instead.)  You can apply it to
your programs, too.

  When we speak of free software, we are referring to freedom, not
price.  Our General Public Licenses are designed to make sure that you
have the freedom to distribute copies of free software (and charge for
this service if you wish), that you receive source code or can get it
if you want it, that you can change the software or use pieces of it
in new free programs; and that you know you can do these things.

  To protect your rights, we need to make restrictions that forbid
anyone to deny you these rights or to ask you to surrender the rights.
These restrictions translate to certain responsibilities for you if you
distribute copies of the software, or if you modify it.

  For example, if you distribute copies of such a program, whether
gratis or for a fee, you must give the recipients all the rights that
you have.  You must make sure that they, too, receive or can get the
source code.  And you must show them these terms so they know their
rights.

  We protect your rights with two steps: (1) copyright the software, and
(2) offer you this license which gives you legal permission to copy,
distribute and/or modify the software.

  Also, for each author's protection and ours, we want to make certain
that everyone understands that there is no warranty for this free
software.  If the software is modified by someone else and passed on, we
want its recipients to know that what they have is not the original, so
that any problems introduced by others will not reflect on the original
authors' reputations.

  Finally, any free program is threatened constantly by software
patents.  We wish to avoid the danger that redistributors of a free
program will individually obtain patent licenses, in effect making the
program proprietary.  To prevent this, we have made it clear that any
patent must be licensed for everyone's free use or not licensed at all.

  The precise terms and conditions for copying, distribution and
modification follow.

   GNU GENERAL PUBLIC LICENSE
   TERMS AND CONDITIONS FOR COPYING, DISTRIBUTION AND MODIFICATION

  0. This License applies to any program or other work which contains
a notice placed by the copyright holder saying it may be distributed
under the terms of this General Public License.  The "Program", below,
refers to any such program or work, and a "work based on the Program"
means either the Program or any derivative work under copyright law:
that is to say, a work containing the Program or a portion of it,
either verbatim or with modifications and/or translated into another
language.  (Hereinafter, translation is included without limitation in
the term "modification".)  Each licensee is addressed as "you".

Activities other than copying, distribution and modification are not
covered by this License; they are outside its scope.  The act of
running the Program is not restricted, and the output from the Program
is covered only if its contents constitute a work based on the
Program (independent of having been made by running the Program).
Whether that is true depends on what the Program does.

  1. You may copy and distribute verbatim copies of the Program's
source code as you receive it, in any medium, provided that you
conspicuously and appropriately publish on each copy an appropriate
copyright notice and disclaimer of warranty; keep intact all the
notices that refer to this License and to the absence of any warranty;
and give any other recipients of the Program a copy of this License
along with the Program.

You may charge a fee for the physical act of transferring a copy, and
you may at your option offer warranty protection in exchange for a fee.

  2. You may modify your copy or copies of the Program or any portion
of it, thus forming a work based on the Program, and copy and
distribute such modifications or work under the terms of Section 1
above, provided that you also meet all of these conditions:

    a) You must cause the modified files to carry prominent notices
    stating that you changed the files and the date of any change.

    b) You must cause any work that you distribute or publish, that in
    whole or in part contains or is derived from the Program or any
    part thereof, to be licensed as a whole at no charge to all third
    parties under the terms of this License.

    c) If the modified program normally reads commands interactively
    when run, you must cause it, when started running for such
    interactive use in the most ordinary way, to print or display an
    announcement including an appropriate copyright notice and a
    notice that there is no warranty (or else, saying that you provide
    a warranty) and that users may redistribute the program under
    these conditions, and telling the user how to view a copy of this
    License.  (Exception: if the Program itself is interactive but
    does not normally print such an announcement, your work based on
    the Program is not required to print an announcement.)

These requirements apply to the modified work as a whole.  If
identifiable sections of that work are not derived from the Program,
and can be reasonably considered independent and separate works in
themselves, then this License, and its terms, do not apply to those
sections when you distribute them as separate works.  But when you
distribute the same sections as part of a whole which is a work based
on the Program, the distribution of the whole must be on the terms of
this License, whose permissions for other licensees extend to the
entire whole, and thus to each and every part regardless of who wrote it.

Thus, it is not the intent of this section to claim rights or contest
your rights to work written entirely by you; rather, the intent is to
exercise the right to control the distribution of derivative or
collective works based on the Program.

In addition, mere aggregation of another work not based on the Program
with the Program (or with a work based on the Program) on a volume of
a storage or distribution medium does not bring the other work under
the scope of this License.

  3. You may copy and distribute the Program (or a work based on it,
under Section 2) in object code or executable form under the terms of
Sections 1 and 2 above provided that you also do one of the following:

    a) Accompany it with the complete corresponding machine-readable
    source code, which must be distributed under the terms of Sections
    1 and 2 above on a medium customarily used for software interchange; or,

    b) Accompany it with a written offer, valid for at least three
    years, to give any third party, for a charge no more than your
    cost of physically performing source distribution, a complete
    machine-readable copy of the corresponding source code, to be
    distributed under the terms of Sections 1 and 2 above on a medium
    customarily used for software interchange; or,

    c) Accompany it with the information you received as to the offer
    to distribute corresponding source code.  (This alternative is
    allowed only for noncommercial distribution and only if you
    received the program in object code or executable form with such
    an offer, in accord with Subsection b above.)

The source code for a work means the preferred form of the work for
making modifications to it.  For an executable work, complete source
code means all the source code for all modules it contains, plus any
associated interface definition files, plus the scripts used to
control compilation and installation of the executable.  However, as a
special exception, the source code distributed need not include
anything that is normally distributed (in either source or binary
form) with the major components (compiler, kernel, and so on) of the
operating system on which the executable runs, unless that component
itself accompanies the executable.

If distribution of executable or object code is made by offering
access to copy from a designated place, then offering equivalent
access to copy the source code from the same place counts as
distribution of the source code, even though third parties are not
compelled to copy the source along with the object code.

  4. You may not copy, modify, sublicense, or distribute the Program
except as expressly provided under this License.  Any attempt
otherwise to copy, modify, sublicense or distribute the Program is
void, and will automatically terminate your rights under this License.
However, parties who have received copies, or rights, from you under
this License will not have their licenses terminated so long as such
parties remain in full compliance.

  5. You are not required to accept this License, since you have not
signed it.  However, nothing else grants you permission to modify or
distribute the Program or its derivative works.  These actions are
prohibited by law if you do not accept this License.  Therefore, by
modifying or distributing the Program (or any work based on the
Program), you indicate your acceptance of this License to do so, and
all its terms and conditions for copying, distributing or modifying
the Program or works based on it.

  6. Each time you redistribute the Program (or any work based on the
Program), the recipient automatically receives a license from the
original licensor to copy, distribute or modify the Program subject to
these terms and conditions.  You may not impose any further
restrictions on the recipients' exercise of the rights granted herein.
You are not responsible for enforcing compliance by third parties to
this License.

  7. If, as a consequence of a court judgment or allegation of patent
infringement or for any other reason (not limited to patent issues),
conditions are imposed on you (whether by court order, agreement or
otherwise) that contradict the conditions of this License, they do not
excuse you from the conditions of this License.  If you cannot
distribute so as to satisfy simultaneously your obligations under this
License and any other pertinent obligations, then as a consequence you
may not distribute the Program at all.  For example, if a patent
license would not permit royalty-free redistribution of the Program by
all those who receive copies directly or indirectly through you, then
the only way you could satisfy both it and this License would be to
refrain entirely from distribution of the Program.

If any portion of this section is held invalid or unenforceable under
any particular circumstance, the balance of the section is intended to
apply and the section as a whole is intended to apply in other
circumstances.

It is not the purpose of this section to induce you to infringe any
patents or other property right claims or to contest validity of any
such claims; this section has the sole purpose of protecting the
integrity of the free software distribution system, which is
implemented by public license practices.  Many people have made
generous contributions to the wide range of software distributed
through that system in reliance on consistent application of that
system; it is up to the author/donor to decide if he or she is willing
to distribute software through any other system and a licensee cannot
impose that choice.

This section is intended to make thoroughly clear what is believed to
be a consequence of the rest of this License.

  8. If the distribution and/or use of the Program is restricted in
certain countries either by patents or by copyrighted interfaces, the
original copyright holder who places the Program under this License
may add an explicit geographical distribution limitation excluding
those countries, so that distribution is permitted only in or among
countries not thus excluded.  In such case, this License incorporates
the limitation as if written in the body of this License.

  9. The Free Software Foundation may publish revised and/or new versions
of the General Public License from time to time.  Such new versions will
be similar in spirit to the present version, but may differ in detail to
address new problems or concerns.

Each version is given a distinguishing version number.  If the Program
specifies a version number of this License which applies to it and "any
later version", you have the option of following the terms and conditions
either of that version or of any later version published by the Free
Software Foundation.  If the Program does not specify a version number of
this License, you may choose any version ever published by the Free Software
Foundation.

  10. If you wish to incorporate parts of the Program into other free
programs whose distribution conditions are different, write to the author
to ask for permission.  For software which is copyrighted by the Free
Software Foundation, write to the Free Software Foundation; we sometimes
make exceptions for this.  Our decision will be guided by the two goals
of preserving the free status of all derivatives of our free software and
of promoting the sharing and reuse of software generally.

       NO WARRANTY

  11. BECAUSE THE PROGRAM IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW.  EXCEPT WHEN
OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
PROVIDE THE PROGRAM "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED
OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE ENTIRE RISK AS
TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU.  SHOULD THE
PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING,
REPAIR OR CORRECTION.

  12. IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES,
INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING
OUT OF THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED TO
LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR
THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS),
EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH
DAMAGES.

       END OF TERMS AND CONDITIONS

*/




/*
---------------------------------------------------------------
ARWEN v1.2.3 Dean Laslett
---------------------------------------------------------------


arwen detects metazoan mitochondrial tRNA genes.

Usage:
arwen -v -s -d -c -l -j -a -q -rn -w -gc -tv
      -seq -br -fasta -fo -o <outfile> <filename>

<filename> is assumed to contain one or more sequences
in FASTA format. Results of the search are printed to
STDOUT. All switches are optional and
case-insensitive.

    -mtmam        Search for Mammalian mitochondrial tRNA
                  genes. -tv switch set.
                  Mammalian mitochondrial genetic code used.
    -mtx          Low scoring tRNA genes are not reported.
    -gc<num>      Use the GenBank transl_table = <num> genetic code.
    -gcstd        Use standard genetic code.
    -gcmet        Use composite Metazoan mitochondrial genetic code.
    -gcvert       Use Vertebrate mitochondrial genetic code.
    -gcinvert     Use Invertebrate mitochondrial genetic code.
    -gcyeast      Use Yeast mitochondrial genetic code.
    -gcprot       Use Mold/Protozoan/Coelenterate mitochondrial genetic code.
    -gcciliate    Use Ciliate genetic code.
    -gcflatworm   Use Echinoderm/Flatworm mitochondrial genetic code
    -gceuplot     Use Euplotid genetic code.
    -gcbact       Use Bacterial/Plant Chloroplast genetic code.
    -gcaltyeast   Use alternative Yeast genetic code.
    -gcascid      Use Ascidian Mitochondrial genetic code.
    -gcaltflat    Use alternative Flatworm Mitochondrial genetic code.
    -gcblep       Use Blepharisma genetic code.
    -gcchloroph   Use Chlorophycean Mitochondrial genetic code.
    -gctrem       Use Trematode Mitochondrial genetic code.
    -gcscen       Use Scenedesmus obliquus Mitochondrial genetic code.
    -gcthraust    Use Thraustochytrium Mitochondrial genetic code.
                  Individual modifications can be appended using\n");
    ,BBB=<aa>     B = A,C,G, or T. <aa> is the three letter\n");
                  code for an amino-acid. More than one modification
                  can be specified. eg -gcvert,aga=Trp,agg=Trp uses
                  the Vertebrate Mitochondrial code and the codons
                  AGA and AGG changed to Tryptophan.          
    -tv           Do not search for mitochondrial TV replacement
                  loop tRNA genes.
    -c7           Search for tRNA genes with 7 base C-loops only.
    -c            Assume that each sequence has a circular
                  topology. Search wraps around each end.
                  Default setting.
    -l            Assume that each sequence has a linear
                  topology. Search does not wrap.
    -d            Double. Search both strands of each
                  sequence. Default setting.
    -s  or -s+    Single. Do not search the complementary
                  (antisense) strand of each sequence.
    -sc or -s-    Single complementary. Do not search the sense
                  strand of each sequence.
    -ps           Lower scoring thresholds to 95% of default levels
    -ps<num>      Change scoring thresholds to <num> percent of default levels
    -rp           Report possible pseudogenes (normalised score < 100).
    -seq          Print out primary sequence.
    -br           Show secondary structure of tRNA gene primary sequence
                  using round brackets.
    -fasta        Print out primary sequence in fasta format.
    -fo           Print out primary sequence in fasta format only 
                  (no secondary structure). 
    -v            Verbose. Prints out information during
                  search to STDERR.
    -a7           Restrict tRNA astem length to a maximum of 7 bases
    -aa           Display message if predicted iso-acceptor species
                  does not match species in sequence name (if present).
    -j            Display 4-base sequence on 3' end of astem
                  regardless of predicted amino-acyl acceptor length.
    -jr           Allow some divergence of 3' amino-acyl acceptor
                  sequence from NCCA.
    -jr4          Allow some divergence of 3' amino-acyl acceptor
                  sequence from NCCA, and display 4 bases.
    -q            Dont print configuration line (which switchs
                  and files were used).
    -rn           Repeat sequence name before summary information.
    -w            Print out in Batch mode.
    -O <outfile>  Print output to <outfile>. If <outfile>
                  already exists, it is overwritten.  By default
                  all output goes to stdout.


*/


#include <stdio.h>
#include <stdlib.h>

#ifndef SEEK_SET
#define SEEK_SET        0
#define SEEK_CUR        1
#define SEEK_END        2
#endif


#define DLIM            '\n'
#define STRLEN          4001
#define STRLENM1        4000
#define KEYLEN          15
#define INACTIVE        2.0e+35
#define IINACTIVE       2000000001L
#define ITHRESHOLD      2000000000L
#define space(c)        (c==' ')||(c=='\t')||(c=='\n')||(c=='\r')
#define sq(pos)         ((pos + d->psmax - 1L) % d->psmax) + 1L
#define itmparam(x,y)   fputc(x,y)

#define noGENE  -1
#define tRNA    0
#define tmRNA   1

#define MAXGCMOD   16
#define MAMMAL_MT  2
#define NGENECODE  24
#define METAZOAN_MT      0
#define STANDARD         1
#define VERTEBRATE_MT    2

#define NAMINOACID 27
#define Phe  0
#define Val  1
#define Leu  2
#define Ile  3
#define Cys  4
#define Gly  5
#define Arg  6
#define Ser  7
#define Ala  8
#define Pro  9
#define Thr  10
#define Tyr  11
#define Asp  12
#define His  13
#define Asn  14
#define Met  15
#define Trp  16
#define Glu  17
#define Gln  18
#define Lys  19
#define Stop 20
#define SeC  21
#define Pyl  22

#define INSERT          -2
#define TERM            -1
#define Adenine         0
#define Cytosine        1
#define Guanine         2
#define Thymine         3
#define AMBIG           4
#define NOBASE          5

#define mtRNAdtthresh   91.5
#define mtRNAtthresh    83.5
#define mtRNAdthresh    85.0
#define PSEUDOGENElevel 0.95




#define RIGHT   0
#define UP      1
#define LEFT    2
#define DOWN    3
#define UPRIGHT 4
#define SLANTDR 5
#define SLANTUR 6
#define SLANTUL 7
#define SLANTDL 8
#define SLANT   5


#define MATX 42
#define MATY 34


#define ASTEM1_EXT       0
#define ASTEM2_EXT       9
#define ASTEM2_EXTD      4                   /* <= ASTEM2_EXT */
#define ASTEM2_EXTE      5                   /* ASTEM2_EXT - ASTEM2_EXTD */
#define MINTSTEM_DIST      (17 + ASTEM2_EXT)
#define MAXTSTEM_DIST      (26 + ASTEM2_EXT)
#define MAXDSTEM_DIST      (9 + ASTEM1_EXT)
#define MINDSTEM_DIST      (8 + ASTEM1_EXT)
#define MININTRONLEN    0
#define MAXINTRONLEN    3000
#define MINCTRNALEN     62
#define MAXCTRNALEN     110
#define MINTRNALEN      (MINCTRNALEN + ASTEM1_EXT + 1)
#define MAXTRNALEN      (MAXCTRNALEN + ASTEM1_EXT + ASTEM2_EXT)
#define MAXETRNALEN     (MAXTRNALEN + MAXINTRONLEN)
#define VARMAX          26 /* 25 */
#define VARMIN          3
#define VARDIFF         23 /* 22 */                   /* VARMAX - VARMIN */
#define MINTPTSDIST     50
#define MAXTPTSDIST     321
#define TPWINDOW        (MAXTPTSDIST - MINTPTSDIST + 1)
#define MINTPDIST       50
#define MAXTPDIST       250
#define TPDISTWINDOW    (MAXTPDIST - MINTPDIST + 1)
#define MINTAGDIST      12
#define MAXTAGDIST      102
#define TAGWINDOW       MAXTAGDIST - MINTAGDIST
#define MINRNACDIST     (MINTPDIST - 5)
#define MAXRNACDIST     (MAXTPDIST - 5)
#define MAXPPINTRONDIST 250
#define TMPTRAILER      145
#define MINPPASDIST     MINTSTEM_DIST
#define MAXPPASDIST     MAXTSTEM_DIST + MAXPPINTRONDIST
#define MINPPTSTPDIST   MINTSTEM_DIST + MINTPDIST
#define MAXPPTSTPDIST   MAXTSTEM_DIST+ASTEM2_EXT+MAXTPDIST+MAXPPINTRONDIST
#define MAXTMRNALEN     (4 + MAXPPASDIST + MAXTPDIST + MAXTAGDIST + TMPTRAILER)
#define TSWEEP          1000
#define WRAP            2*MAXETRNALEN
#define NPTAG           33

/*
NOTE: If MAXPPINTRONDIST is increased, then validity of MAXTMRNALEN
and MAXETRNALEN must be ensured. WRAP = 2*MAXETRNALEN determines the length
of wseq, which contains the wrap around for circular sequences. This
must remain equal to or more than 2*MAXTMRNALEN and TSWEEP.
*/

#define BASE   0
#define FSTEM  1
#define BSTEM  2

#define NOID   0
#define DLOOP  1
#define DSTEM  2
#define CLOOP  3
#define VAR    4

#define NA   MAXINTRONLEN
#define ND   100
#define NT   150  /* 150, 500 */
#define NH   2000
#define NTH  3000
#define NC   5000 
#define NS   5
#define NTAG 474  /* 367 */
#define LSEQ 20000
#define ATBOND 2.5
#define mtNA 1500
#define mtND 150 
#define mtNTH 750 
#define mtNTM  3
#define mtGCBOND    0.0
#define mtATBOND    -0.5
#define mtGTBOND    -1.2
#define mtTTBOND    -2.9
#define mtGGBOND    -3.0
#define mtGABOND    -3.0
#define mtNOBOND    -3.0
#define mtBONDSTAB  1.5
#define mtABONDSTAB 2.0
#define mtTSTTSTAB   -2.5
#define mtTERMSTAB  0.01
#define mtSENDSTAB  0.01
#define mtNSTAB     0.1
#define mt3MMSTAB   1.0
#define mtGCPENALTY 0.8
#define mtGCPENALTYD 2.0
#define mt_DRLmaxlength 16
#define mt_TVRLmaxlength 18
#define mtNCLM 3


typedef struct { char filename[80];
                 FILE *f;
                 char seqname[STRLEN];
                 double gc;
                 long ps;
                 long psmax;
                 long seqstart;
                 long nextseq;
                 int ns; } data_set;


typedef struct { char name[80];
                 int seq[MAXTRNALEN+1];
                 int eseq[MAXETRNALEN+1];
                 int *ps;
                 int nbase;
                 int comp;
                 long start;
                 long stop;
                 int astem1;
                 int astem2;
                 int aatail;
                 int spacer1;
                 int spacer2;
                 int dstem;
                 int dloop;
                 int cstem;
                 int cloop;
                 int intron;
                 int nintron;
                 int anticodon;
                 int var;
                 int varbp;
                 int tstem;
                 int tloop;
                 int genetype;
                 double energy;
                 int asst;
                 int tps;
                 int tpe;   } gene;

typedef struct { int *pos;
                 int stem;
                 int loop;
                 double energy; } trna_loop;

typedef struct { int *pos;
                 int stem;
                 int loop;
                 unsigned int bondtype;
                 double energy;
                 double stem_energy; } mt_trna_loop;

typedef struct { int *pos;
                 int *looppos;
                 int *end;
                 int stem;
                 int loop;
                 int arm;
                 int anticodon;
                 unsigned int bondtype;
                 double energy;
                 double stem_energy; } mt_trna_cloop;

typedef struct { int *pos;
                 int stem;
                 int loop;
		 int *end;
                 unsigned int bondtype;
                 double energy;
                 double stem_energy; } mt_trna_tloop;


typedef struct { int *pos;
                 int *end;
                 int stem;
                 int loop;
                 double energy; } trna_dloop;


typedef struct { int *pos1;
                 int *pos2;
                 int stem;
                 double energy; } trna_astem;

typedef struct { int *pos1;
                 int *pos2;
                 int stem;
                 unsigned int bondtype;
                 double energy; } mt_trna_astem;


typedef struct { char arg[50];
                 FILE *f;
                 int batch;
                 int repeatsn;
                 int trna;
                 int mtrna;
                 int tvloop;
                 int cloop7;
                 int geneticcode;
                 int ngcmod;
                 int gcmod[MAXGCMOD];
                 int gcfix;
                 int discrim;
                 int extastem;
                 int tarm;
                 int tagthresh;
                 int tarmlength;
                 int showconfig;
                 int libflag;
                 int verbose;
                 int linear;
                 int both;
                 int reportpseudogenes;
                 int energydisp;
                 int trnadomaindisp;
                 int seqdisp;
                 int aataildisp;
                 int aataildiv;
                 int sp1max;
		         int mtxdetect;
		         int matchacceptor;
                 double trnathresh;
                 double tmrnathresh;
                 double tarmthresh;
                 double tcthresh;
                 double tthresh;
                 double dthresh;
                 double tathresh;
                 double tmathresh;
                 double tmcthresh;
                 double tmcathresh;
                 int tmstrict;
                 unsigned int tmrthresh;
                 unsigned int itathresh;
                 int maxintronlen;
                 int minintronlen;
                 int ioverlay;
		         int ifixedpos;
                 double mttthresh;
                 double mtdthresh;
                 double mtdtthresh;
                 double mttarmthresh;
                 double mtdarmthresh;
                 int ngene[NS];
		         int iamismatch;
                 int loffset;
                 int roffset;
                 long start;
                 int comp;
                 int space;
                 int trna_struct[200];
               } csw;


/* Basepair matching matrices */


  int bp[6][6] = { { 0,0,0,1,1,0 },
                   { 0,0,1,0,1,0 },
                   { 0,1,0,1,1,0 },
                   { 1,0,1,0,1,0 },
                   { 1,1,1,1,1,0 },
                   { 0,0,0,0,0,0 } };

  int wbp[6][6] =
   { { 0,0,0,2,2,0 },
     { 0,0,2,0,2,0 },
     { 0,2,0,1,2,0 },
     { 2,0,1,0,2,0 },
     { 2,2,2,2,2,0 },
     { 0,0,0,0,0,0 } };

  int wcbp[6][6] = { { 0,0,0,1,1,0 },
                     { 0,0,1,0,1,0 },
                     { 0,1,0,0,1,0 },
                     { 1,0,0,0,1,0 },
                     { 1,1,1,1,1,0 },
                     { 0,0,0,0,0,0 } };


  int gc[6][6] = { { 0,0,0,0,0,0 },
                   { 0,0,1,0,1,0 },
                   { 0,1,0,0,1,0 },
                   { 0,0,0,0,0,0 },
                   { 0,1,1,0,1,0 },
                   { 0,0,0,0,0,0 } };

  int gt[6][6] = { { 0,0,0,0,0,0 },
                   { 0,0,0,0,0,0 },
                   { 0,0,0,1,1,0 },
                   { 0,0,1,0,1,0 },
                   { 0,0,1,1,1,0 },
                   { 0,0,0,0,0,0 } };

  int at[6][6] = { { 0,0,0,1,1,0 },
                   { 0,0,0,0,0,0 },
                   { 0,0,0,0,0,0 },
                   { 1,0,0,0,0,0 },
                   { 1,0,0,0,1,0 },
                   { 0,0,0,0,0,0 } };

  int tt[6][6] = { { 0,0,0,0,0,0 },
                   { 0,0,0,0,0,0 },
                   { 0,0,0,0,0,0 },
                   { 0,0,0,1,1,0 },
                   { 0,0,0,1,1,0 },
                   { 0,0,0,0,0,0 } };

  int stemterm[6][6] = { { 0,0,1,0,1,0 },
                         { 0,0,0,0,0,0 },
                         { 1,0,0,0,1,0 },
                         { 0,0,0,1,1,0 },
                         { 1,0,1,1,1,0 },
                         { 0,0,0,0,0,0 } };

  int aastemterm[6][6] =
   { { 1,0,1,0,1,0 },
     { 0,0,0,0,0,0 },
     { 1,0,0,0,1,0 },
     { 0,0,0,1,1,0 },
     { 1,0,1,1,1,0 },
     { 0,0,0,0,0,0 } };

  int ggstemterm[6][6] =
   { { 0,0,1,0,1,0 },
     { 0,0,0,0,0,0 },
     { 1,0,1,0,1,0 },
     { 0,0,0,1,1,0 },
     { 1,0,1,1,1,0 },
     { 0,0,0,0,0,0 } };

  int assymst[6][6] = { { 0,0,0,0,0,0 },
                        { 0,0,0,0,0,0 },
                        { 1,0,0,0,1,0 },
                        { 0,0,0,1,1,0 },
                        { 1,0,0,1,1,0 },
                        { 0,0,0,0,0,0 } };

  int assymat[6][6] = { { 0,0,0,1,1,0 },
                        { 0,0,0,0,0,0 },
                        { 0,0,0,0,0,0 },
                        { 0,0,0,0,0,0 },
                        { 0,0,0,1,1,0 },
                        { 0,0,0,0,0,0 } };


  int stackbp[6][6] = { { 0,0,0,1,1,0 },
                        { 0,0,1,0,1,0 },
                        { 0,1,0,1,1,0 },
                        { 1,0,1,1,1,0 },
                        { 1,1,1,1,1,0 },
                        { 0,0,0,0,0,0 } };

  int ggstackbp[6][6] =
   { { 0,0,0,1,1,0 },
     { 0,0,1,0,1,0 },
     { 0,1,1,1,1,0 },
     { 1,0,1,1,1,0 },
     { 1,1,1,1,1,0 },
     { 0,0,0,0,0,0 } };


  int ggbp[6][6] =
   { { 0,0,0,1,1,0 },
     { 0,0,1,0,1,0 },
     { 0,1,1,1,1,0 },
     { 1,0,1,0,1,0 },
     { 1,1,1,1,1,0 },
     { 0,0,0,0,0,0 } };

  int gabp[6][6] =
   { { 0,0,1,1,1,0 },
     { 0,0,1,0,1,0 },
     { 1,1,0,1,1,0 },
     { 1,0,1,0,1,0 },
     { 1,1,1,1,1,0 },
     { 0,0,0,0,0,0 } };

  int assymagbp[6][6] =
   { { 0,0,1,1,1,0 },
     { 0,0,1,0,1,0 },
     { 0,1,0,1,1,0 },
     { 1,0,1,0,1,0 },
     { 1,1,1,1,1,0 },
     { 0,0,0,0,0,0 } };

  int stembp[6][6] =
   { { 0,0,1,1,1,0 },
     { 0,0,1,0,1,0 },
     { 1,1,0,1,1,0 },
     { 1,0,1,1,1,0 },
     { 1,1,1,1,1,0 },
     { 0,0,0,0,0,0 } };

  int ggstembp[6][6] =
   { { 0,0,1,1,1,0 },
     { 0,0,1,0,1,0 },
     { 1,1,1,1,1,0 },
     { 1,0,1,1,1,0 },
     { 1,1,1,1,1,0 },
     { 0,0,0,0,0,0 } };

  int gastembp[6][6] =
   { { 1,0,1,1,1,0 },
     { 0,0,1,0,1,0 },
     { 1,1,1,1,1,0 },
     { 1,0,1,1,1,0 },
     { 1,1,1,1,1,0 },
     { 0,0,0,0,0,0 } };


  int vbp[6][6] =
   { { 0,0,1,4,4,0 },
     { 0,0,4,0,4,0 },
     { 1,4,0,2,4,0 },
     { 4,0,2,0,4,0 },
     { 4,4,4,4,4,0 },
     { 0,0,0,0,0,0 } };


  int tandemid[mtNTM][4] =
   { { 3,2,2,3 },
     { 2,3,3,2 },
     { 3,3,3,3 } };

  double tandem_em[mtNTM] = { -0.5,-0.5,2.0 };


  double send_em[6][6] =
   { { 0.0,0.0,0.0,0.0,0.0,0.0 },
     { 0.0,0.0,0.5*mtSENDSTAB,0.0,0.5*mtSENDSTAB,0.0 },
     { 0.0,0.5*mtSENDSTAB,0.0,mtSENDSTAB,mtSENDSTAB,0.0 },
     { 0.0,0.0,mtSENDSTAB,0.0,mtSENDSTAB,0.0 },
     { 0.0,0.5*mtSENDSTAB,mtSENDSTAB,mtSENDSTAB,mtSENDSTAB,0.0 },
     { 0.0,0.0,0.0,0.0,0.0,0.0 } };


  double ssend_em[6][6] =
   { { 0.0,0.0,0.0,0.0,0.0,0.0 },
     { 0.0,0.0,mtSENDSTAB,0.0,mtSENDSTAB,0.0 },
     { 0.0,mtSENDSTAB,0.0,mtSENDSTAB,mtSENDSTAB,0.0 },
     { 0.0,0.0,mtSENDSTAB,0.0,mtSENDSTAB,0.0 },
     { 0.0,mtSENDSTAB,mtSENDSTAB,mtSENDSTAB,mtSENDSTAB,0.0 },
     { 0.0,0.0,0.0,0.0,0.0,0.0 } };



  int neighbour_map[6][6] =
   { { 0,0,1,0,1,0 },
     { 0,0,0,0,0,0 },
     { 1,0,0,0,1,0 },
     { 0,0,0,1,1,0 },
     { 1,0,1,1,1,0 },
     { 0,0,0,0,0,0 } };

  double neighbour_em[2][6][6] = {
   { { 0.0,0.0,0.0,0.0,0.0,0.0 },
     { 0.0,0.0,0.0,0.0,0.0,0.0 },
     { 0.0,0.0,0.0,0.0,0.0,0.0 },
     { 0.0,0.0,0.0,0.0,0.0,0.0 },
     { 0.0,0.0,0.0,0.0,0.0,0.0 },
     { 0.0,0.0,0.0,0.0,0.0,0.0 } },

   { { 0.0,0.0,0.0,0.0,0.0,0.0 },
     { 0.0,0.0,mtNSTAB,0.0,mtNSTAB,0.0 },
     { 0.0,mtNSTAB,0.0,0.0,mtNSTAB,0.0 },
     { 0.0,0.0,0.0,0.0,0.0,0.0 },
     { 0.0,mtNSTAB,mtNSTAB,0.0,mtNSTAB,0.0 },
     { 0.0,0.0,0.0,0.0,0.0,0.0 } } };

  unsigned int btmap[6][6] =
   { { 0x10000,0x10000,0x1000,0x10,0x00000,0x10000 },
     { 0x10000,0x10000,0x1,0x10000,0x00000,0x10000 },
     { 0x1000,0x1,0x10000,0x100,0x00000,0x10000 },
     { 0x10,0x10000,0x100,0x1000,0x00000,0x10000 },
     { 0x00000,0x00000,0x00000,0x00000,0x00000,0x10000 },
     { 0x10000,0x10000,0x10000,0x10000,0x10000,0x10000 } };



  int mt_discrim[3][64][6] =
   /* metazoan mt */
   {{{ 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },

     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 0,0,0,0,0,0 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },

     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },

     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 0,0,0,0,0,0 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 }},
  /* standard */
    {{ 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },

     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },

     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },

     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 }},
   /* mammal mt */
    {{ 1,0,0,0,1,1 },  { 1,0,0,0,1,1 },  { 1,0,0,0,1,1 },  { 1,0,0,0,1,1 },
     { 0,0,0,1,1,1 },  { 1,0,0,0,1,1 },  { 1,0,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,0,1,0,1,1 },  { 1,0,0,0,1,1 },  { 1,1,1,1,1,1 },
     { 1,0,0,0,1,1 },  { 1,1,1,1,1,1 },  { 0,1,0,0,1,1 },  { 0,0,1,0,1,1 },

     { 1,0,0,0,1,1 },  { 1,0,0,0,1,1 },  { 1,0,0,0,1,1 },  { 1,0,1,0,1,1 },
     { 0,0,1,0,1,1 },  { 1,0,0,0,1,1 },  { 1,0,1,1,1,1 },  { 1,0,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,0,1,0,1,1 },  { 1,0,0,0,1,1 },  { 1,1,1,1,1,1 },
     { 0,0,0,0,0,0 },  { 1,0,1,1,1,1 },  { 0,0,1,0,1,1 },  { 1,1,1,1,1,1 },

     { 1,0,0,0,1,1 },  { 1,0,0,0,1,1 },  { 1,0,0,0,1,1 },  { 1,0,1,0,1,1 },
     { 0,1,0,1,1,1 },  { 1,0,0,0,1,1 },  { 1,0,1,1,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,1,1,1,1 },  { 1,0,1,0,1,1 },  { 1,0,0,0,1,1 },  { 1,1,1,1,1,1 },
     { 1,1,0,0,1,1 },  { 1,1,1,1,1,1 },  { 0,1,0,0,1,1 },  { 0,0,1,0,1,1 },

     { 1,1,0,1,1,1 },  { 1,0,0,0,1,1 },  { 1,0,1,1,1,1 },  { 1,0,0,0,1,1 },
     { 1,0,1,0,1,1 },  { 1,0,0,0,1,1 },  { 1,1,1,1,1,1 },  { 0,0,0,0,0,0 },
     { 1,1,1,1,1,1 },  { 1,0,1,0,1,1 },  { 1,0,1,0,1,1 },  { 1,1,1,1,1,1 },
     { 1,0,0,0,1,1 },  { 1,1,1,1,1,1 },  { 1,0,1,0,1,1 },  { 1,1,1,1,1,1 }}};

/* GENETIC CODES (INDEXED BY ANTICODON) */

  char aapolarity[NAMINOACID+1] = "NNNNPNPPNNPNPPPNNPPP***????";
  char aaletter[NAMINOACID+1] = "FVLICGRSAPTYDHNMWEQK***????";
  char aaname[NAMINOACID][20] =
   { "Phe","Val","Leu","Ile","Cys",
     "Gly","Arg","Ser","Ala","Pro",
     "Thr","Tyr","Asp","His","Asn",
     "Met","Trp","Glu","Gln","Lys",
     "Stop",
     "seC",
     "Pyl",
     "(Arg|Stop|Ser|Gly)",
     "(Ile|Met)",
     "(Stop|Trp)",
     "(Lys|Asn)" };

  char ambig_aaname[4] = "???";

  int aamap[NGENECODE][64] = {
   /* composite metazoan mt */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,23,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,24,
     25,Gly,Arg,23,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,26 },
   /* standard */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     SeC,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
  /* vertebrate mt */
  { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Stop,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Stop,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
   /* yeast mt */
   { Phe,Val,Thr,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Thr,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Thr,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Thr,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
   /* mold, protozoan, and coelenterate mt */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
   /* invertebrate mt */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
   /* ciliate */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Gln,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     SeC,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Gln,Glu,Gln,Lys },
   /* deleted -> standard */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     SeC,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
   /* deleted -> standard */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     SeC,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
   /* echinoderm and flatworm mt */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     Trp,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Asn },
   /* euplotid */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     Cys,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
   /* bacterial and plant chloroplast */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Ser,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     SeC,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
   /* alternate yeast */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Ser,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     SeC,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
   /* ascidian mt */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Gly,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Gly,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
   /* alternate flatworm mt */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     Trp,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Glu,Gln,Asn },
   /* blepharisma */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Gln,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     SeC,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
   /* chlorophycean mt */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Leu,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     SeC,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
   /* deleted -> standard */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     SeC,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
   /* deleted -> standard */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     SeC,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
   /* deleted -> standard */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     SeC,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
   /* deleted -> standard */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     SeC,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
   /* trematode mt */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
   /* scenedesmus obliquus mt*/
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Leu,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Leu,Ile,
     SeC,Gly,Arg,Arg,
     Stop,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys },
   /* thraustochytrium mt */
   { Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Leu,Val,Ser,Met,
     Trp,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Pyl,Glu,Gln,Lys,
     Phe,Val,Leu,Ile,
     Cys,Gly,Arg,Ser,
     Ser,Ala,Pro,Thr,
     Tyr,Asp,His,Asn,
     Stop,Val,Leu,Ile,
     SeC,Gly,Arg,Arg,
     Ser,Ala,Pro,Thr,
     Stop,Glu,Gln,Lys } };


/* POINTERS TO DETECTED GENES */

  gene *ts;

/* TOOLS */

char upcasec(char c)
{ return((c >= 'a')?c-32:c); }


int length(char *s)
{ int i=0;
  while (*s++) i++;
  return(i); }

char *strpos(char *s, char *k)
{ char c,d;
  int i;
  d = *k;
  while (c = *s)
   { if (c == d)
       { i = 0;
         do if (!k[++i]) return(s); while (s[i] == k[i]); }
     s++; }
  return(NULL); }


char *softstrpos(char *s, char *k)
{ char c,d;
  int i;
  d = upcasec(*k);
  while (c = *s)
   { if (upcasec(c) == d)
       { i = 0;
         do if (!k[++i]) return(s);
         while (upcasec(s[i]) == upcasec(k[i])); }
     s++; }
  return(NULL); }



char *dconvert(char *s, double *r)
{ static char zero='0',nine='9';
  int shift,expshift,sgn,expsgn,exponent;
  char c,limit;
  double result;
  shift = 0;
  expshift = 0;
  sgn = 1;
  expsgn = 1;
  limit = 0;
  exponent = 0;
  result = 0.0;
  if ((c = *s) == '-')
   { sgn = -1;
     c = *++s; }
  else if (c == '+') c= *++s;
  if (c >= zero)
   if (c <= nine)
    { result = (double)(c - zero);
      while ((c = *++s) >= zero)
       { if (c > nine) break;
         if (++limit < 15) result = result*10.0 + (double)(c - zero); }}
  if (c == '.')
   while ((c = *++s) >= zero)
    { if (c > nine) break;
      if (++limit < 15)
       { result = result*10.0 + (double)(c - zero);
         shift++; }}
  if ((c == 'E')||(c == 'e')||(c == 'D')||(c == 'd'))
    { if ((c = *++s) == '-')
       { expsgn = -1;
         c = *++s; }
      else
       if (c == '+') c = *++s;
      if (c >= zero)
       if (c <= nine)
         { exponent = c - zero;
           while ((c = *++s) >= zero)
            { if (c > nine) break;
              exponent = exponent*10 + c - zero;
              if (++expshift > 3) break; }}}
  result *= (double)sgn;
  exponent = exponent*expsgn - shift;
  if (exponent >= 0)
    while (exponent--) result *= 10.0;
  else
    while (exponent++) result /= 10.0;
  (*r) *= 0.01*result;
  return(s); }

char *iconvert(char *s, int *r)
{ static char zero='0',nine='9';
  int sgn;
  int result;
  char c;
  sgn = 1;
  result = 0;
  if ((c = *s) == '-')
   { sgn = -1;
     c = *++s; }
  else if (c == '+') c= *++s;
  if (c >= zero)
   if (c <= nine)
    { result = (int)(c - zero);
      while ((c = *++s) >= zero)
       { if (c > nine) break;
         result = result*10 + (int)(c - zero); }}
  *r = result * sgn;
  return(s); }



char *copy(char *from, char *to)
{ while (*to++ = *from++);
  return(--to);  }


char *copy3cr(char *from, char *to, int n)
{ while (*to = *from++)
   { if (*to == DLIM)
      { *to = '\0';
        break; }
     if (--n <= 0)
      { *++to = '\0';
        break; }
     to++; }
  return(to); }

/* LIBRARY */

long get_seqname(data_set *d)
{ int i,ic;
  char c,*s;
  FILE *f;
  f = d->f;
  fseek(f,d->seqstart,SEEK_SET);
  do { if ((ic = getc(f)) == EOF) return(-1L);
       c = (char)ic; }
  while (space(c));
  if (c != '>')
   { *d->seqname = '\0';
     fseek(f,d->seqstart,SEEK_SET);
     return(d->seqstart); }
  else
   { if (!fgets(d->seqname,STRLENM1,f)) return(-1L);
     s = d->seqname;
     i = 0;
     while ((c = *s) != '\0')
      { if (c == '\n') break;
        if (c == '\r') break;
        if (++i >= STRLEN) break;
        s++; }
     *s = '\0';
     return(ftell(f)); }}


int move_forward(data_set *d)
{ int ic;
  static int map[256] =
  { -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,
    -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,
    -4,-4,-4,-4,-4,-4,NOBASE,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,
    -4,-4,-2,-4,-4,Adenine,AMBIG,Cytosine,AMBIG,-4,-4,Guanine,AMBIG,
    -4,-4,AMBIG,-4,AMBIG,AMBIG,-4,-4,-4,
    AMBIG,AMBIG,Thymine,Thymine,AMBIG,AMBIG,-4,
    AMBIG,-4,-4,-4,-4,INSERT,NOBASE,-4,Adenine,AMBIG,Cytosine,AMBIG,
    -4,-4,Guanine,AMBIG,-4,-4,AMBIG,-4,AMBIG,AMBIG,-4,-4,-4,
    AMBIG,AMBIG,Thymine,Thymine,AMBIG,AMBIG,-4,
    AMBIG,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,
    -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,
    -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,
    -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,
    -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,
    -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,
    -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4 };
  if (d->ps >= d->psmax)
   if (d->psmax > 0L)
    { fseek(d->f,d->seqstart,SEEK_SET);
      d->ps = 0L; }
  NL:
  if ((ic = getc(d->f)) == EOF) goto FAIL;
  ic = map[ic];
  if (ic >= Adenine)
   { d->ps++;
     return(ic); }
  if (ic == -2)
   { d->nextseq = ftell(d->f) - 1L;
     return(TERM); }
  goto NL;
  FAIL:
  d->nextseq = -1L;
  if (d->psmax > 0L)
   { d->ps = d->psmax;
     return(NOBASE); }
  else return(TERM); }


int seq_init(data_set *d)
{ long ngc;
  int ic;
  if ((d->seqstart = get_seqname(d)) < 0L) return(0);
  d->ps = 0L;
  d->psmax = -1L;
  ngc = 0L;
  while ((ic = move_forward(d)) >= Adenine)
   if (ic >= Cytosine)
    if (ic <= Guanine)
     ngc++;
  if ((d->psmax = d->ps) <= 0L) return(0);
  d->gc = (double)ngc/(double)d->psmax;
  fseek(d->f,d->seqstart,SEEK_SET);
  d->ps = 0L;
  return(1); }



char cbase(int c)
{ static char base[6] = "acgt..";
  if (c < Adenine) return('#');
  if (c > NOBASE) return((char)c);
  return(base[c]); }


char cpbase(int c)
{ static char base[6] = "ACGT..";
  if (c < Adenine) return('#');
  if (c > NOBASE) return((char)c);
  return(base[c]); }


char *aa(int *anticodon, csw *sw)
{ int p1,p2,p3;
  if ((p1 = *anticodon) >= AMBIG) return(ambig_aaname);
  if ((p2 = anticodon[1]) >= AMBIG) return(ambig_aaname);
  if ((p3 = anticodon[2]) >= AMBIG) return(ambig_aaname);
  return(aaname[aamap[sw->geneticcode][(p1<<4) + (p2<<2) + p3]]); }


double gc_content(gene *t)
{ int *s,*se;
  double ngc;
  static double score[6] = { 0.0,1.0,1.0,0.0,0.0,0.0 };
  ngc = 0.0;
  if ((t->nintron > 0) && (t->asst == 0))
   { s = t->eseq;
     se = s + t->intron;
     while (s < se) ngc += score[*s++];
     s = se + t->nintron;
     se = t->eseq + t->nbase + t->nintron;
     while (s < se) ngc += score[*s++]; }
  else
   { s = t->seq;
     se = s + t->nbase;
     while (s < se) ngc += score[*s++]; }
  return(ngc/(double)t->nbase); }


void write_seq(FILE *f, int *seq, int newline)
{ int i,c;
  i = 0;
  while ((c = *seq++) >= Adenine)
   { fputc(cbase(c),f);
     if (newline)
      if (++i >= 50)
       { fputc('\n',f);
         i = 0; }}
  if (i > 0) fputc('\n',f); }


int find_var_hairpin(gene *t)
{ int e,stem,vstem,loop,*sn,*sen,*pos1,*pos2,*sb,*se,*sc,*sd,*sf,*s;
  unsigned int c,cn,m;
  static unsigned int A[6] = { 0,0,0x100,0x400,0,0 };
  static unsigned int C[6] = { 0,0,0x400,0,0,0 };
  static unsigned int G[6] = { 0x100,0x400,0,0x200,0,0 };
  static unsigned int T[6] = { 0x400,0,0x200,0,0,0 };
  static unsigned int te[6] = { 0,0,0,0,0,0 };
  if (t->genetype != tRNA) return(0);
  if (t->var < 13) return(0);
  e = 0;
  sb = t->seq + t->astem1 + t->spacer1 + 2*t->dstem + t->dloop + 
       t->spacer2 + 2*t->cstem + t->cloop + t->nintron;
  sc = sb + 3;  /* 4 */
  se = sb + t->var - 2;  /* 3 */
  sf = se - 2;
  te[0] = A[*se];
  te[1] = C[*se];
  te[2] = G[*se];
  te[3] = T[*se];
  while (--se > sf)
   { te[0] = (te[0] >> 4) | A[*se];
     te[1] = (te[1] >> 4) | C[*se];
     te[2] = (te[2] >> 4) | G[*se];
     te[3] = (te[3] >> 4) | T[*se]; }
  while (se >= sc)
   { te[0] = ((te[0] >> 4) | A[*se]);
     te[1] = ((te[1] >> 4) | C[*se]);
     te[2] = ((te[2] >> 4) | G[*se]);
     te[3] = ((te[3] >> 4) | T[*se]);
     s = se - 5;
     sd = se - 7;
     m = te[*s];
     while (--s > sd) m = (m >> 4) + te[*s];
     while (s >= sb)
       {  m = (m >> 4) + te[*s];
          c = m & 0xf;
          if (c >= 9)
           { stem = 3;
             loop = (int)(se - s) - 3;
             sen = se;
             sn = s + 2;
             while (loop >= 6)
              { if ((cn = vbp[sen[-1]][sn[1]]) <= 0) break;
                c += cn;
                stem++;
                loop -= 2;
                sen--;
                sn++; }
             if (c > e)
              { e = c;
                pos1 = s;
                pos2 = sen;
                vstem = stem; }}
          s--; }
      se--; }
  if (e > 0)
   return((((int)(pos1 - sb)) << 10) + (((int)(pos2 - sb)) << 5) + vstem); 
  else
   return(0); }    

void write_to_library(FILE *f, gene *t, csw *sw)
{ int *s;
  static char trnatype[2][6] = { "tRNA","mtRNA" };
  s = t->seq + t->anticodon;
  fprintf(f,">%s",t->name);
  if (!softstrpos(t->name,"RNA"))
   switch (t->genetype)
    { case tRNA:
      default:
       t->varbp = find_var_hairpin(t);
       if (t->tstem == 0) fprintf(f," TV-loop");
       else if (t->dstem == 0) fprintf(f," D-loop");
       switch(t->cloop)
        { case 6:
           fprintf(f," %s-???(%c%c)",trnatype[sw->mtrna],
                     cbase(*s),cbase(*(s+1)));
           break;
          case 8:
           fprintf(f," %s-???(%c%c%c%c)",trnatype[sw->mtrna],
                     cbase(*s),cbase(s[1]),cbase(s[2]),cbase(s[3]));
           break;
          case 7:
          default:
           fprintf(f," %s-%s(%c%c%c)",trnatype[sw->mtrna],
                     aa(s,sw),cbase(*s),cbase(*(s+1)),cbase(*(s+2)));
           break; }
      break; }
  if (strpos(t->name,"bases)"))
   fprintf(f,"\n");
  else
   fprintf(f," (%d bases)\n",t->nbase);
  fprintf(f,"sequence =\n");
  write_seq(f,t->seq,1);
  if (*t->eseq >= Adenine)
   { fprintf(f,"extended sequence =\n");
     write_seq(f,t->eseq,1); }
  fprintf(f,"nbase = %d\n",t->nbase);
  fprintf(f,"sense = %d\n",t->comp);
  fprintf(f,"start = %ld\n",t->start);
  fprintf(f,"stop = %ld\n",t->stop);
  fprintf(f,"astem1 = %d\n",t->astem1);
  fprintf(f,"astem2 = %d\n",t->astem2);
  fprintf(f,"atail = %d\n",t->aatail);
  fprintf(f,"spacer1 = %d\n",t->spacer1);
  fprintf(f,"spacer2 = %d\n",t->spacer2);
  fprintf(f,"dstem = %d\n",t->dstem);
  fprintf(f,"dloop = %d\n",t->dloop);
  fprintf(f,"cstem = %d\n",t->cstem);
  fprintf(f,"cloop = %d\n",t->cloop);
  fprintf(f,"anticodon = %d\n",t->anticodon);
  fprintf(f,"nintron = %d\n",t->nintron);
  fprintf(f,"intron = %d\n",t->intron);
  fprintf(f,"var = %d\n",t->var);
  fprintf(f,"varbp = %d,%d,%d\n",((t->varbp >> 10)&0x1f),
            ((t->varbp >> 5)&0x1f),(t->varbp&0x1f));
  fprintf(f,"tstem = %d\n",t->tstem);
  fprintf(f,"tloop = %d\n",t->tloop);
  fprintf(f,"gc = %g\n",gc_content(t)); }



void init_trna(FILE *f, csw *sw)
{ int c,*s;
  s = sw->trna_struct;
  while ((c = *s++) != TERM) itmparam(cbase(c),f); }




int *make_tv(int *seq, char matrix[][MATY],
             int *x, int *y, int orient, int tv)
{ int i,px,py,stem;
  static int ux[4] = { 1,0,-1,0 };
  static int uy[4] = { 0,1,0,-1 };
  static int vx[4] = { 0,-1,0,1 };
  static int vy[4] = { 1,0,-1,0 };
  static int loopu[26][26] =
  { { 0 }, { 0 }, { 0 }, { 0 }, { 0 }, { 0 },
    { 1,1,1,0,0,-1,-1 },
    { 1,1,1,1,0,-1,-1,-1 },
    { 1,1,1,1,0,0,-1,-1,-1 },
    { 1,1,1,1,1,0,-1,-1,-1,-1 },
    { 1,1,1,1,1,0,0,-1,-1,-1,-1 },
    { 1,1,1,1,1,1,0,-1,-1,-1,-1,-1 },
    { 1,1,1,1,1,1,0,0,-1,-1,-1,-1,-1 },
    { 1,1,1,1,1,1,0,0,0,-1,-1,-1,-1,-1 },
    { 1,1,1,1,1,1,1,0,0,-1,-1,-1,-1,-1,-1 },
    { 1,1,1,1,1,1,1,0,0,0,-1,-1,-1,-1,-1,-1 },
    { 1,1,1,1,1,1,1,1,0,0,-1,-1,-1,-1,-1,-1,-1 },
    { 1,1,1,1,1,1,1,1,0,0,0,-1,-1,-1,-1,-1,-1,-1 },
    { 1,1,1,1,1,1,1,1,1,0,0,-1,-1,-1,-1,-1,-1,-1,-1 },
    { 1,1,1,1,1,1,1,1,1,0,0,0,-1,-1,-1,-1,-1,-1,-1,-1 },
    { 1,1,1,1,1,1,1,1,1,1,0,0,-1,-1,-1,-1,-1,-1,-1,-1,-1 },
    { 1,1,1,1,1,1,1,1,1,1,0,0,0,-1,-1,-1,-1,-1,-1,-1,-1,-1 },
    { 1,1,1,1,1,1,1,1,1,1,1,0,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1 },
    { 1,1,1,1,1,1,1,1,1,1,1,0,0,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1 },
    { 1,1,1,1,1,1,1,1,1,1,1,1,0,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1 },
    { 1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1 } };
  static int loopv[26][26] =
  { { 0 }, { 0 }, { 0 }, { 0 }, { 0 }, { 0 },
    { 0,1,1,1,1,1,1 },
    { -1,1,1,1,1,1,1,1 },
    { -1,1,1,1,1,1,1,1,0 },
    { -1,0,1,1,1,1,1,1,1,0 },
    { -1,0,1,1,1,1,1,1,1,0,0 },
    { -1,0,0,1,1,1,1,1,1,1,0,0 },
    { -1,0,0,1,1,1,1,1,1,1,0,0,0 },
    { -1,0,0,1,1,1,1,1,1,1,0,0,0,0 },
    { -1,0,0,1,1,1,1,1,1,1,0,0,0,0,0 },
    { -1,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0 },
    { -1,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0 },
    { -1,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0 },
    { -1,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0 },
    { -1,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0 },
    { -1,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0 },
    { -1,0,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0 },
    { -1,0,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0 },
    { -1,0,0,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0 },
    { -1,0,0,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0 },
    { -1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0 } };
  px = *x;
  py = *y;
  stem = 0;
  if (tv < 6)
   { px += ux[orient];
     py += uy[orient];
     i = 0;
     while (i < tv)
     { px += vx[orient];
       py += vy[orient];
       matrix[px][py] = cbase(*seq++);
       i++; }
     py += (6-i)*vy[orient];
     goto FN; }
  if (tv > 25)
   { if (tv % 2)
      stem = (tv - 25)/2;
     else
      stem = (tv - 24)/2;
     tv = tv - 2*stem; }
  i = 0;
  while (i < stem)
   { px += ux[orient];
     py += uy[orient];
     matrix[px][py] = cbase(*seq++);
     i++; }
  i = 0;
  while (i < tv)
  { px += ux[orient]*loopu[tv][i] + vx[orient]*loopv[tv][i];
    py += uy[orient]*loopu[tv][i] + vy[orient]*loopv[tv][i];
    matrix[px][py] = cbase(*seq++);
    i++; }
  px += ux[orient]*loopu[tv][i] + vx[orient]*loopv[tv][i];
  py += uy[orient]*loopu[tv][i] + vy[orient]*loopv[tv][i];
  i = 0;
  while (i < stem)
   { matrix[px][py] = cbase(*seq++);
     px -= (ux[orient]);
     py -= (uy[orient]);
     i++; }
  FN:
  *x = px;
  *y = py;
  return(seq); }


int base_match(char b1, char b2)
{ int i,s;
  static char base1[11] = "acgtgtagtg";
  static char base2[11] = "tgcatggatg";
  static int score[11] = { 2,2,2,2,1,1,3,3,3,3 };
  s = 0;
  for (i = 0; i < 10; i++)
   if (b1 == base1[i])
    if (b2 == base2[i])
     { s = score[i];
       break; }
  return(s); }


int *make_clover(int *seq, int b, int e, int stemlength,
                  char matrix[][MATY], int *x, int *y, int orient)
{ int i,px,py,pxb,pyb,pxe,pye,l,xlg,xlgd,ylgh,ylg;
  int *s,*se;
  static int ux[9] = { 1,0,-1,0,0,1,1,-1,-1 };
  static int uy[9] = { 0,1,0,-1,1,-1,1,1,-1 };
  static int vx[9] = { 0,-1,0,1,1,1,1,-1,-1 };
  static int vy[9] = { 1,0,-1,0,0,0,0,0,0 };
  static int loopu[18][18] =
  { { -1 }, { 0,-1 }, { 0,0,-1 }, { 0,1,-1,-1 }, { 0,1,0,-1,-1 },
    { 0,1,0,0,-1,-1 }, { 0,1,1,0,-1,-1,-1 }, { 0,1,1,0,0,-1,-1,-1 },
    { 0,1,1,1,0,-1,-1,-1,-1 }, { 0,1,1,1,0,0,-1,-1,-1,-1 },
    { 0,1,1,1,0,0,0,-1,-1,-1,-1 },
    { 0,1,1,1,1,0,0,-1,-1,-1,-1,-1 },
    { 0,1,1,1,1,0,0,0,-1,-1,-1,-1,-1 },
    { 0,1,1,1,1,1,0,0,-1,-1,-1,-1,-1,-1 },
    { 0,1,1,1,1,1,0,0,0,-1,-1,-1,-1,-1,-1 },
    { 0,1,1,1,1,1,1,0,0,-1,-1,-1,-1,-1,-1,-1 },
    { 0,1,1,1,1,1,1,0,0,0,-1,-1,-1,-1,-1,-1,-1 },
    { 0,1,1,1,1,1,1,0,0,0,0,-1,-1,-1,-1,-1,-1,-1 } };
  static int loopv[18][18] =
  { { 2 }, { 1,1 },  { 0,1,1 }, { -1,2,2,-1 }, { -1,1,1,2,-1 },
    { -1,1,1,1,1,-1 }, { -1,0,1,1,1,1,-1 }, { -1,0,1,1,1,1,0,-1 },
    { -1,0,1,1,1,1,0,0,-1 }, { -1,0,0,1,1,1,1,0,0,-1 },
    { -1,0,0,0,1,1,1,1,0,0,-1 },
    { -1,0,0,0,1,1,1,1,0,0,0,-1 },
    { -1,0,0,0,1,1,1,1,0,0,0,0,-1 },
    { -1,0,0,0,0,1,1,1,1,0,0,0,0,-1 },
    { -1,0,0,0,0,1,1,1,1,0,0,0,0,0,-1 },
    { -1,0,0,0,0,0,1,1,1,1,0,0,0,0,0,-1 },
    { -1,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,-1 },
    { -1,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,-1 } };
  static int dloopu[18][18] =
  { { -1 }, { 0,-1 }, { 0,0,-1 }, { 0,1,-1,-1 }, { 0,1,0,-1,-1 },
    { 0,1,0,0,-1,-1 }, { 0,1,1,0,-1,-1,-1 }, { 0,1,1,0,0,-1,-1,-1 },
    { 0,1,1,0,0,0,-1,-1,-1 }, { 0,1,1,1,0,0,-1,-1,-1,-1 },
    { 0,1,1,1,0,0,0,-1,-1,-1,-1 },
    { 0,1,1,1,1,0,0,-1,-1,-1,-1,-1 },
    { 0,1,1,1,1,0,0,0,-1,-1,-1,-1,-1 },
    { 0,1,1,1,1,1,0,0,-1,-1,-1,-1,-1,-1 },
    { 0,1,1,1,1,1,0,0,0,-1,-1,-1,-1,-1,-1 },
    { 0,1,1,1,1,1,1,0,0,-1,-1,-1,-1,-1,-1,-1 },
    { 0,1,1,1,1,1,1,0,0,0,-1,-1,-1,-1,-1,-1,-1 },
    { 0,1,1,1,1,1,1,0,0,0,0,-1,-1,-1,-1,-1,-1,-1 } };
  static int dloopv[18][18] =
  { { 2 }, { 1,1 },  { 0,1,1 }, { -1,2,2,-1 }, { -1,1,1,2,-1 },
    { -1,1,1,1,1,-1 }, { -1,0,1,1,1,1,-1 }, { -1,0,1,1,1,1,0,-1 },
    { -1,0,1,1,1,1,1,-1,-1 }, { -1,0,0,1,1,1,1,0,0,-1 },
    { -1,0,0,0,1,1,1,1,0,0,-1 },
    { -1,0,0,0,1,1,1,1,0,0,0,-1 },
    { -1,0,0,0,1,1,1,1,0,0,0,0,-1 },
    { -1,0,0,0,0,1,1,1,1,0,0,0,0,-1 },
    { -1,0,0,0,0,1,1,1,1,0,0,0,0,0,-1 },
    { -1,0,0,0,0,0,1,1,1,1,0,0,0,0,0,-1 },
    { -1,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,-1 },
    { -1,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,-1 } };
  static char bond1[5] = " +!:";
  static char bond2[5] = " +-.";
  px = *x;
  py = *y;
  s = seq + b;
  se = s + stemlength;
  while (s  < se)
  { matrix[px][py] = cbase(*s++);
    px += ux[orient];
    py += uy[orient]; }
  l = e - b - 2*stemlength;
  if (l < 0) l = 0;
  if (l < 18)
   { i = 0;
     if (orient == DOWN)
      { while (i < l)
         { px += ux[orient]*dloopu[l][i] + vx[orient]*dloopv[l][i];
           py += uy[orient]*dloopu[l][i] + vy[orient]*dloopv[l][i];
           matrix[px][py] = cbase(*s++);
           i++; }
         px += ux[orient]*dloopu[l][i] + vx[orient]*dloopv[l][i];
         py += uy[orient]*dloopu[l][i] + vy[orient]*dloopv[l][i]; }
     else
      { while (i < l)
         { px += ux[orient]*loopu[l][i] + vx[orient]*loopv[l][i];
           py += uy[orient]*loopu[l][i] + vy[orient]*loopv[l][i];
           matrix[px][py] = cbase(*s++);
           i++; }
         px += ux[orient]*loopu[l][i] + vx[orient]*loopv[l][i];
         py += uy[orient]*loopu[l][i] + vy[orient]*loopv[l][i]; }}
  else 
   { ylgh = ((l >> 2) - 2) >> 1;
     ylg = (ylgh << 1) + 2;
     xlgd = l - ylg - 2*ylgh + 1;
     xlg = (xlgd + 1) >> 1;
     pxb = px - ylgh*vx[orient];
     if ((pxb < 0) || (pxb >= MATX)) goto NOLOOP;
     pyb = py - ylgh*vy[orient];
     if ((pyb < 0) || (pyb >= MATY)) goto NOLOOP;
     pxe = px + xlg*ux[orient] + (ylg - ylgh + 1)*vx[orient];
     if ((pxe < 0) || (pxe >= MATX)) goto NOLOOP;
     pye = py + xlg*uy[orient] + (ylg - ylgh + 1)*vy[orient];
     if ((pye < 0) || (pye >= MATY)) goto NOLOOP;  
     for (i = 0; i < ylgh; i++)
      { px -= vx[orient];
        py -= vy[orient];
        matrix[px][py] = cbase(*s++); }
     for (i = 0; i < xlg; i++)
      { px += ux[orient];
        py += uy[orient];
        matrix[px][py] = cbase(*s++); }
     for (i = 1; i < ylg; i++)
      { px += vx[orient];
        py += vy[orient];
        matrix[px][py] = cbase(*s++); }
     px += vx[orient];
     py += vy[orient];
     if (!(xlgd & 1)) matrix[px][py] = cbase(*s++);
     for (i = 0; i < xlg; i++)
      { px -= ux[orient];
        py -= uy[orient];
        matrix[px][py] = cbase(*s++); }
     for (i = 1; i < ylgh; i++)
      { px -= vx[orient];
        py -= vy[orient];
        matrix[px][py] = cbase(*s++); }
     px -= (ux[orient] + vx[orient]);
     py -= (uy[orient] + vy[orient]); }
  goto STEMBOND;
  NOLOOP:
  px += ux[orient]*loopu[0][0] + vx[orient]*loopv[0][0];
  py += uy[orient]*loopu[0][0] + vy[orient]*loopv[0][0];
  STEMBOND:
  se = seq + e;
  s = se - stemlength;
  while (s  < se)
  { matrix[px][py] = cbase(*s++);
    i = base_match(matrix[px][py],
                    matrix[px - 2*vx[orient]][py - 2*vy[orient]]);
    switch(orient)
     { case  RIGHT:
       case  LEFT:  matrix[px - vx[orient]][py - vy[orient]] = bond1[i];
                    break;
       case  SLANTDR:
       case  SLANTUR:
       case  SLANTUL:
       case  SLANTDL:
       case  UPRIGHT:
       case  UP:
       case  DOWN:  matrix[px - vx[orient]][py - vy[orient]] = bond2[i];
                    break; }
    px -= ux[orient];
    py -= uy[orient]; }
  *x = px;
  *y = py;
  return(se); }



int *make_dv(int *seq, char matrix[][MATY], int dloop,
                  int orient, int *xp, int *yp)
{ int i,x,y;
  static int ux[5] = { 1,0,-1,0,0 };
  static int uy[5] = { 0,1,0,-1,1 };
  static int vx[5] = { 0,-1,0,1,1 };
  static int vy[5] = { 1,0,-1,0,0 };
  static int loopu[22][22] =
  { { -1 }, { -1,0 },
    { -1,-1,1 },
    { -1,-1,0,1 },
    { -1,-1,0,0,1 },
    { -1,-1,-1,0,1,1 },
    { -1,-1,-1,0,0,1,1 },
    { -1,-1,-1,-1,0,1,1,1 },
    { -1,-1,-1,-1,0,0,1,1,1 },
    { -1,-1,-1,-1,-1,0,1,1,1,1 },
    { -1,-1,-1,-1,-1,0,0,1,1,1,1 },
    { -1,-1,-1,-1,-1,-1,0,1,1,1,1,1 },
    { -1,-1,-1,-1,-1,-1,0,0,1,1,1,1,1 },
    { -1,-1,-1,-1,-1,-1,-1,0,1,1,1,1,1,1 },
    { -1,-1,-1,-1,-1,-1,0,-1,0,1,1,1,1,1,1 },
    { -1,-1,-1,-1,-1,-1,-1,-1,0,1,1,1,1,1,1,1 },
    { -1,-1,-1,-1,-1,-1,-1,-1,0,0,1,1,1,1,1,1,1 },
    { -1,-1,-1,-1,-1,-1,-1,-1,-1,0,1,1,1,1,1,1,1,1 },
    { -1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,1,1,1,1,1,1,1,1 },
    { -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,1,1,1,1,1,1,1,1,1 },
    { -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,1,1,1,1,1,1,1,1,1 },
    { -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,1,1,1,1,1,1,1,1,1,1 } };
  static int loopv[22][22] =
  { { -6 }, { -3,-3 },
    { -2,-2,-2 },
    { -2,-1,-1,-2 },
    { -1,-1,-2,-1,-1 },
    { -1,-1,-1,-1,-1,-1 },
    { 0,-1,-1,-1,-1,-1,-1 },
    { 0,-1,0,-1,-1,-1,-1,-1 },
    { 0,-1,0,-1,-1,-1,0,-1,-1 },
    { 0,0,-1,0,-1,-1,-1,0,-1,-1 },
    { 0,0,-1,0,-1,-1,-1,0,-1,0,-1 },
    { 0,0,0,-1,0,-1,-1,-1,0,-1,0,-1 },
    { 0,0,0,-1,0,-1,-1,-1,0,-1,0,0,-1 },
    { 0,0,0,0,-1,0,-1,-1,-1,0,-1,0,0,-1 },
    { 0,0,0,0,-1,0,-1,-1,-1,0,-1,0,0,0,-1 },
    { 0,0,0,0,0,-1,0,-1,-1,-1,0,-1,0,0,0,-1 },
    { 0,0,0,0,0,-1,0,-1,-1,-1,-1,0,0,0,0,0,-1 },
    { 0,0,0,0,0,0,-1,0,-1,-1,-1,0,-1,0,0,0,0,-1 },
    { 0,0,0,0,0,0,-1,0,-1,-1,-1,-1,0,0,0,0,0,0,-1 },
    { 0,0,0,0,0,0,0,-1,0,-1,-1,-1,0,-1,0,0,0,0,0,-1 },
    { 0,0,0,0,0,0,0,-1,0,-1,-1,-1,-1,0,0,0,0,0,0,0,-1 },
    { 0,0,0,0,0,0,0,0,-1,0,-1,-1,-1,0,-1,0,0,0,0,0,0,-1 } };
  x = *xp;
  y = *yp;
  if ((dloop < 2) || (dloop > 21))
  { x--;
    y-= 6;
    seq += dloop;
    goto FN; }
  i = 0;
  while (i < dloop)
  { x += ux[orient]*loopu[dloop][i] + vx[orient]*loopv[dloop][i];
    y += uy[orient]*loopu[dloop][i] + vy[orient]*loopv[dloop][i];
    matrix[x][y] = cbase(*seq++);
    i++; }
  x += ux[orient]*loopu[dloop][i] + vx[orient]*loopv[dloop][i];
  y += uy[orient]*loopu[dloop][i] + vy[orient]*loopv[dloop][i];
  FN:
  *xp = x;
  *yp = y;
  return(seq); }



int *make_var(int *seq, char matrix[][MATY],
               int *x, int *y, int orient, int var, int varbp)
{ int i,b,e,p,px,py,pxf,pyf,l,stem;
  static int ux[4] = { 1,0,-1,0 };
  static int uy[4] = { 0,1,0,-1 };
  static int vx[4] = { 0,-1,0,1 };
  static int vy[4] = { 1,0,-1,0 };
  static int preu[4][5][4] =
   { { {},{1},{1},{1},{1} },
     { {},{1,1},{1,1},{1,1},{1,1} },
     { {},{0,1,1},{0,1,1},{1,1,1},{1,1,1} },
     { {},{0,0,1,1},{0,0,1,1},{0,1,1,1},{0,1,1,1} } };
  static int prev[4][5][4] =
   { { {},{0},{0},{0},{-1} },
     { {},{0,1},{0,0},{-1,0},{-1,0} },
     { {},{0,1,0},{0,0,0},{-1,0,0},{-1,0,-1} },
     { {},{0,1,0,0},{0,1,0,-1},{0,0,0,-1},{0,0,-1,-1} } };
  static int postu[4][5][4] =
   { { {},{0},{1,-1},{0,-1,0},{1,-1,-1,0} },
     { {},{0},{0,-1},{0,-1,0},{0,-1,-1,0} },
     { {},{0},{0,-1},{0,-1,-1},{1,-1,-1,-1} },
     { {},{0},{0,-1},{0,-1,-1},{1,-1,-1,-1} } };
  static int postv[4][5][4] =
   { { {},{0},{0,1},{0,1,1},{0,1,1,1} },
     { {},{0},{0,1},{0,1,1},{0,1,1,1} },
     { {},{0},{0,1},{0,1,1},{0,1,1,1} },
     { {},{0},{0,1},{0,1,1},{0,1,1,1} } };
  static int loopu[10][10] =
  { { 2 }, { 1,1 }, { 1,1,0 }, { 1,0,1,0 }, { 1,1,0,0,0 },
    { 1,1,1,-1,-1,1 }, { 1,1,1,0,0,-1,0 }, { 1,1,1,1,0,-1,-1,0 },
        { 1,1,1,1,0,0,-1,-1,0 }, { 1,1,1,1,1,0,-1,-1,-1,0 } };
  static int loopv[10][10] =
  { { 3 }, { 1,2 },  { 0,1,2 }, { 0,1,1,1 }, { -1,1,1,1,1 },
    { -1,0,1,1,1,1 }, { -1,-1,1,1,1,1,1 }, { -1,-1,0,1,1,1,1,1 },
        { -1,-1,-1,1,1,1,1,1,1 }, { -1,-1,-1,0,1,1,1,1,1,1 } };
  px = *x;
  py = *y;
  if (var < 0) var = 0;
  if (var > 30) var = 30;
  if (varbp > 0)
   { b = (varbp >> 10) & 0x1f;
     if (b > 3) goto NBP;
     stem = varbp & 0x1f;
     e = stem + ((varbp >> 5) & 0x1f);
     p = var - e;
     if (p < 1) goto NBP;  /* 2 */
     if (p > 4) goto NBP;
     pxf = px + 2*ux[orient] + 3*vx[orient];
     pyf = py + 2*uy[orient] + 3*vy[orient];
     i = 0;
     while (i < b)
      { px += ux[orient]*preu[b][p][i] + vx[orient]*prev[b][p][i];
        py += uy[orient]*preu[b][p][i] + vy[orient]*prev[b][p][i];
        matrix[px][py] = cbase(*seq++);
        i++; }
     px += ux[orient]*preu[b][p][b] + vx[orient]*prev[b][p][b];
     py += uy[orient]*preu[b][p][b] + vy[orient]*prev[b][p][b];
     seq = make_clover(seq,0,e-b,stem,matrix,&px,&py,orient+SLANT);
     i = 0;
     while (i < p)
      { px += ux[orient]*postu[b][p][i] + vx[orient]*postv[b][p][i];
        py += uy[orient]*postu[b][p][i] + vy[orient]*postv[b][p][i];
        matrix[px][py] = cbase(*seq++);
        i++; }
     *x = pxf;
     *y = pyf;
     goto FIN;  }
  NBP:   
  if (var > 9)
   { if (var % 2) stem = (var - 7)/2;
     else stem = (var - 6)/2; }
  else stem = 0;
  l = var - 2*stem;
  i = 0;
  while (i < stem)
   { px += ux[orient] - vx[orient];
     py += uy[orient] - vy[orient];
     matrix[px][py] = cbase(*seq++);
     i++; }
  i = 0;
  while (i < l)
   { px += ux[orient]*loopu[l][i] + vx[orient]*loopv[l][i];
     py += uy[orient]*loopu[l][i] + vy[orient]*loopv[l][i];
     matrix[px][py] = cbase(*seq++);
     i++; }
  px += ux[orient]*loopu[l][i] + vx[orient]*loopv[l][i];
  py += uy[orient]*loopu[l][i] + vy[orient]*loopv[l][i];
  i = 0;
  while (i < stem)
   { matrix[px][py] = cbase(*seq++);
     px -= (ux[orient] - vx[orient]);
     py -= (uy[orient] - vy[orient]);
     i++; }
  *x = px;
  *y = py;
  FIN:
  return(seq); }



void remove_inserts(int *s1, int *s2)
{ int flag,c;
  flag = 0;
  while ((c = *s1++) != TERM)
   { if (c == INSERT)
      { flag = 1 - flag;
        continue; }
     if (flag) continue;
     *s2++ = c; }
  *s2 = TERM; }



void build_trna(gene *t, char matrix[][MATY], int x, int y, csw *sw)
{ int i,j,e,c,*seq;
  int rseq[150];
  static char bond2[5] = " +-.";
  t->varbp = find_var_hairpin(t);
  remove_inserts(t->seq,rseq);
  seq = rseq;
  i = 0;
  while (i < t->astem1)
  { matrix[x][y] = cbase(*seq++);
    y--;
    i++; }
  if (t->spacer1 > 0)
   { x--;
     if (t->spacer1 >= 3) matrix[x][y+1] = cbase(*seq++);
     matrix[x][y] = cbase(*seq++);
     y--;
     x--;
     if (t->spacer1 >= 2) matrix[x][y] = cbase(*seq++);
     if ((t->spacer2 < 2) || (t->spacer1 > 1))
      { x--;
        y--; }}
  if (t->dstem > 0)
   { e = 2*t->dstem + t->dloop;
     seq = make_clover(seq,0,e,t->dstem,matrix,&x,&y,LEFT);
     if (t->spacer2 > 1) x--;
     y--;
     if (t->spacer2 > 0) matrix[x][y] = cbase(*seq++);
     y--;
     if (t->spacer2 > 1)
      { if (t->spacer1 > 1) x++;
        matrix[x][y] = cbase(*seq++);
        if (t->spacer1 < 2) y--; }
     x++; }
  else
   seq = make_dv(seq,matrix,t->dloop,RIGHT,&x,&y);
  e = 2*t->cstem + t->cloop;
  seq = make_clover(seq,0,e,t->cstem,matrix,&x,&y,DOWN);
  if (t->tstem > 0)
   { seq = make_var(seq,matrix,&x,&y,RIGHT,t->var,t->varbp);
     e = 2*t->tstem + t->tloop;
     seq = make_clover(seq,0,e,t->tstem,matrix,&x,&y,RIGHT);
     y++; }
  else
   seq = make_tv(seq,matrix,&x,&y,RIGHT,t->tloop);
  e = t->astem2;
  i = 0;
  while (i < e)
  { if ((c = *seq++) < Adenine) break;
    matrix[x][y] = cbase(c);
    j = base_match(matrix[x][y],matrix[x - 2][y]);
    matrix[x - 1][y] = bond2[j];
    y++;
    i++; }
  i = 0;
  e = (sw->aataildisp)?ASTEM2_EXTD:t->aatail;
  j = (e < 2)?e:2;
  while (i < j)
  { if ((c = *seq++) < Adenine) break;
    matrix[x][y] = cbase(c);
    x++;
    y++;
    i++; }
  e -= j;
  i = 0;
  while (i < e)
  { if ((c = *seq++) < Adenine) break;
    matrix[x][y] = cbase(c);
    x++;
    i++; }}






void init_matrix(char matrix[][MATY])
{ int i,j;
  for (i =0; i < MATY; i++)
   for (j = 0; j < MATX; j++) matrix[j][i] = ' '; }


void disp_matrix(FILE *f, char matrix[][MATY], int ylines)
{ int i,j,k;
  i = ylines;
  while (--i >= 0)
   { k = MATX;
     while (--k > 0) if (matrix[k][i] != ' ') break;
     for (j = 0; j <= k; j++) fputc(matrix[j][i],f);
     fputc('\n',f); }
  fputc('\n',f); }


void xcopy(char m[][MATY], int x, int y, char *s, int l)
{ int i;
  char c;
  i = 0;
  while (i < l)
   { if (x >= MATX) break;
     if (!(c = *s++)) break;
     m[x++][y] = c;
     i++; }}





void sense_switch(int *seq1, int *seq2, int lseq)
{ int i,b;
  int *sseq,*cseq;
  sseq = seq1;
  cseq = seq2 + lseq;
  while (--cseq >= seq2)
   { b = *sseq++;
     if (b >= Adenine)
      { if (b <= Thymine)
         *cseq = Thymine - b;
        else
         { if (b <= NOBASE) *cseq = b;
           else *cseq = NOBASE; }}
     else *cseq = NOBASE; }}



char *position(char *s, gene *t, csw *sw)
{ long start;
  start = t->start;
  if (sw->linear) if (start <= 0) start--;
  if (t->comp)
   sprintf(s,"c[%ld,%ld]",start,t->stop);
  else
   sprintf(s,"[%ld,%ld]",start,t->stop);
  return(s); }


void location(char *s, gene *t, csw *sw, char *m)
{ char sp[80];
  sprintf(s,"%s %s",m,position(sp,t,sw)); }

void disp_location(gene *t, csw *sw, char *m)
{ char sp[80];
  fprintf(sw->f,"%s %s\n",m,position(sp,t,sw)); }

char *name(gene *t, char *si, int proc, csw *sw)
{ int *s;
  char *sb,*st;
  static char trnatype[2][6] = { "tRNA","mtRNA" };
  switch (t->genetype)
   { case tRNA:
              s = (proc?t->seq:t->ps) + t->anticodon;
              st = trnatype[sw->mtrna];
              sb = si;
              if (t->dstem == 0)
               { sprintf(sb,"D-loop ");
                 sb += 7; }
              if (t->tstem == 0)
               { sprintf(sb,"TV-loop ");
                 sb += 8; }
              if (t->cloop == 8)
               sprintf(sb,"%s-?(%s|%s)(%c%c%c%c)",st,
			  aa(s,sw),aa(s+1,sw),
                          cbase(*s),cbase(s[1]),cbase(s[2]),
                          cbase(s[3]));
              else if (t->cloop == 6)
               sprintf(sb,"%s-?(%s|%s)(%c%c)",st,
			  aa(s-1,sw),aa(s,sw),
                          cbase(*s),cbase(s[1]));
              else
               sprintf(sb,"%s-%s(%c%c%c)",st,
                          aa(s,sw),cbase(*s),cbase(s[1]),cbase(s[2]));
              break;
     default: *si = '\0';
              break; }
  return(si); }


void disp_fasta_seq(FILE *f, gene *t, int ns, int n, int nsp, int c, csw *sw)
{ int i,*s,*se;
  char genename[100],genepos[100];
  if (t->nintron > 0)
   { s = t->eseq;
     se = s + t->nbase + t->nintron; }
  else
   { s = t->seq;
     se = s + t->nbase; }
  name(t,genename,1,sw);
  position(genepos,t,sw);
  if (nsp > 0)
   { if (ns > 0) fprintf(f,">%d-%d%s%s\n",ns,n,genename,genepos);
     else fprintf(f,">%s%s\n",genename,genepos); }
  else
   { if (ns > 0) fprintf(f,">%d-%d %s %s\n",ns,n,genename,genepos);
     else fprintf(f,">%s %s\n",genename,genepos); }
  i = 0;
  while (s < se)
   { if (c) fputc(cpbase(*s++),f);
     else fputc(cbase(*s++),f);
     if (++i >= 50)
      { fputc('\n',f);
        i = 0; }}
  if (i > 0) fputc('\n',f); }


void disp_seq(FILE *f, gene *t, csw *sw)
{ int i,j,k,varbp,stem,ab,ae,hl,*s,*se,*sl,*sb,*sr;
  char genename[100];
  static int bplb[2] = { '.','(' };
  static int bprb[2] = { '.',')' };
  if (t->nintron > 0)
   { s = t->eseq;
     se = s + t->nbase + t->nintron; }
  else
   { s = t->seq;
     se = s + t->nbase; }
  if (sw->seqdisp >= 3)
   { if (!sw->batch) fputc('\n',f);
     if (sw->seqdisp == 3) disp_fasta_seq(f,t,0,0,0,0,sw);
     else disp_fasta_seq(f,t,0,0,0,1,sw); }
  else
   { if (!sw->batch)
      { name(t,genename,1,sw);
        fprintf(f,"\nPrimary sequence for %s\n",genename); }
     if ((sw->seqdisp == 2) && (t->genetype == tRNA))
     { sl = s;
       while (sl < se) fputc(cbase(*sl++),f);
       fputc('\n',f);
       sl = s;
       sr = se - t->aatail - 1;
       for (i = 0; i < t->astem1; i++,sl++,sr--) fputc(bplb[bp[*sl][*sr]],f);
       for (i = 0; i < t->spacer1; i++) fputc(' ',f);
       sl += t->spacer1;
       sb = sl + t->dstem - 1;
       sr = sb + t->dstem + t->dloop;      
       for (i = 0; i < t->dstem; i++,sl++,sr--) fputc(bplb[bp[*sl][*sr]],f);
       for (i = 0; i < t->dloop; i++) fputc('d',f);
       sl += t->dloop;
       for (i = 0; i < t->dstem; i++,sl++,sb--) fputc(bprb[bp[*sl][*sb]],f);
       for (i = 0; i < t->spacer2; i++) fputc(' ',f);
       sl += t->spacer2;
       sb = sl + t->cstem - 1;
       sr = sb + t->cstem + t->cloop + t->nintron;
       for (i = 0; i < t->cstem; i++,sl++,sr--) fputc(bplb[bp[*sl][*sr]],f);
       hl = t->astem1 + t->spacer1 + 2*t->dstem + t->dloop + 
              t->spacer2 + t->cstem;
       if (t->nintron > 0)
        { j = t->intron - hl;
          ab = t->anticodon - hl;
          ae = ab + t->cloop - 5;
          for (i = 0; i < j; i++) 
           if (i <= ae)
            if (i >= ab) fputc('A',f);
            else fputc(' ',f);
           else fputc(' ',f);
          for (i = 0; i < t->nintron; i++) fputc('i',f);
          for (i = j; i < t->cloop; i++) 
           if (i <= ae)
            if (i >= ab) fputc('A',f);
            else fputc(' ',f);
           else fputc(' ',f); }
       else   
        { j = t->cloop - 4;
          ab = t->anticodon - hl;
          ae = t->cloop - ab - j;
          for (i = 0; i < ab; i++) fputc(' ',f);
          for (i = 0; i < j; i++) fputc('A',f);
          for (i = 0; i < ae; i++) fputc(' ',f); }
       sl += (t->cloop + t->nintron);
       for (i = 0; i < t->cstem; i++,sl++,sb--) fputc(bprb[bp[*sl][*sb]],f);
       varbp = find_var_hairpin(t);
       if (varbp > 0)
        { j = (varbp >> 10);
          k = (varbp >> 5) & 0x1f;
          stem = (varbp & 0x1f);
          sr = sl + k + stem - 1;
          sl += j;
          sb = sl + stem - 1;
          for (i = 0; i < j; i++) fputc(' ',f);
          for (i = 0; i < stem; i++,sl++,sr--) fputc(bplb[bp[*sl][*sr]],f);
          for (i = j+stem; i < k; i++,sl++) fputc('v',f);
          for (i = 0; i < stem; i++,sl++,sb--) fputc(bprb[bp[*sl][*sb]],f);
          for (i = k+stem; i < t->var; i++,sl++) fputc(' ',f); }
       else
        { for (i = 0; i < t->var; i++) fputc(' ',f);
          sl += t->var; }
       sb = sl + t->tstem - 1;
       sr = sb + t->tstem + t->tloop;
       for (i = 0; i < t->tstem; i++,sl++,sr--) fputc(bplb[bp[*sl][*sr]],f);
       for (i = 0; i < t->tloop; i++) fputc('t',f);
       sl += t->tloop;
       for (i = 0; i < t->tstem; i++,sl++,sb--) fputc(bprb[bp[*sl][*sb]],f);
       sb = s + t->astem1 - 1;
       for (i = 0; i < t->astem2; i++,sl++,sb--) fputc(bprb[bp[*sl][*sb]],f);
       fputc('\n',f); }
     else
     { if (!sw->batch)
         fprintf(f,"1   .   10    .   20    .   30    .   40    .   50\n");
       i = 0;
       while (s < se)
        { fputc(cbase(*s++),f);
          if (++i >= 50)
           { fputc('\n',f);
             i = 0; }}
       if (i > 0) fputc('\n',f); }}
  if (!sw->batch) 
   { fputc('\n',f);
     fputc('\n',f); }}





double nenergy(gene *t, csw *sw)
{ double eref;
  eref = 100.0;
  if (t->genetype == tRNA)
   if (sw->mtrna)
    { if (t->dstem == 0) eref = mtRNAtthresh;
      else 
       if (t->tstem == 0) eref = mtRNAdthresh;
       else eref = mtRNAdtthresh; }
  return(100.0*t->energy/eref); }



void disp_gene(gene *t, char m[][MATY], csw *sw)
{ double gc;
  char stat[80];
  switch(t->genetype)
   { case tRNA:
             build_trna(t,m,13,27,sw);
             name(t,stat,1,sw);
             xcopy(m,4,3,stat,length(stat));
             break; }
  location(stat,t,sw,"Sequence");
  xcopy(m,4,1,stat,length(stat));
  gc = gc_content(t);
  sprintf(stat,"%d bases, %%GC = %2.1f",t->nbase,100.0*gc);
  xcopy(m,4,2,stat,length(stat));
  if (sw->reportpseudogenes)
   if (t->energy < 100.0) 
    xcopy(m,4,4,"Possible Pseudogene",19);
  if (sw->energydisp)
   { sprintf(stat,"Score = %g\n",t->energy);
     xcopy(m,4,0,stat,length(stat)); }}

void disp_batch_trna(FILE *f, gene *t, csw *sw)
{ int ls,ps,*s,anticodon;
  char pos[50],species[50];
  static char type[2][6] = { "tRNA","mtRNA" };
  static char asterisk[2] = { ' ','*'};
  anticodon = 1 + t->anticodon;
  if (t->nintron > 0)
   if (t->intron <= t->anticodon)
    anticodon += t->nintron;
  s = t->seq + t->anticodon;
  ps = sw->reportpseudogenes?((t->energy < 100.0)?1:0):0;
  switch(t->cloop)
   { case 6:
     case 8:
	  sprintf(species,"%s-???%c",type[sw->mtrna],asterisk[ps]);
	  break;
     case 7:
	 default:
          sprintf(species,"%s-%s%c",type[sw->mtrna],aa(s,sw),asterisk[ps]);
	  break; }
  position(pos,t,sw);
  ls = length(species);
  if (ls <= 10) fprintf(f,"%-10s%28s",species,pos);
  else if (ls <= 17) fprintf(f,"%-17s%21s",species,pos);
  else fprintf(f,"%-25s%13s",species,pos);
  if (sw->energydisp)
   fprintf(f,"\t%5.1f",t->energy);
  fprintf(f,"\t%-4d",anticodon);
  switch(t->cloop)
   { case 6:
          fprintf(f,"\t(%c%c) ",cbase(*s),cbase(s[1]));
	  break;
     case 8:
          fprintf(f,"\t(%c%c%c%c) ",
	          cbase(*s),cbase(s[1]),cbase(s[2]),cbase(s[3]));
	  break;
     case 7:
     default:
          fprintf(f,"\t(%c%c%c)",cbase(*s),cbase(s[1]),cbase(s[2]));
	  break; }
  if (t->nintron > 0)
   fprintf(f,"i(%d,%d)",t->intron+1,t->nintron);
  fputc('\n',f);
  if (sw->seqdisp) disp_seq(f,t,sw); }



double vloop_stability(int *sb, int var, int *varbp)
{ int e,stem,vstem,loop,*sn,*sen,*pos1,*pos2,*se,*sc,*sd,*sf,*s;
  unsigned int c,cn,m;
  static unsigned int A[6] = { 0,0,0x100,0x400,0,0 };
  static unsigned int C[6] = { 0,0,0x400,0,0,0 };
  static unsigned int G[6] = { 0x100,0x400,0,0x200,0,0 };
  static unsigned int T[6] = { 0x400,0,0x200,0,0,0 };
  static unsigned int te[6] = { 0,0,0,0,0,0 };
  e = 0;
  sc = sb + 3;   
  se = sb + var - 2; 
  sf = se - 2;
  te[0] = A[*se];
  te[1] = C[*se];
  te[2] = G[*se];
  te[3] = T[*se];
  while (--se > sf)
   { te[0] = (te[0] >> 4) | A[*se];
     te[1] = (te[1] >> 4) | C[*se];
     te[2] = (te[2] >> 4) | G[*se];
     te[3] = (te[3] >> 4) | T[*se]; }
  while (se >= sc)
   { te[0] = ((te[0] >> 4) | A[*se]);
     te[1] = ((te[1] >> 4) | C[*se]);
     te[2] = ((te[2] >> 4) | G[*se]);
     te[3] = ((te[3] >> 4) | T[*se]);
     s = se - 5;
     sd = se - 7;
     m = te[*s];
     while (--s > sd) m = (m >> 4) + te[*s];
     while (s >= sb)
       {  m = (m >> 4) + te[*s];
          c = m & 0xf;
          if (c >= 9)
           { stem = 3;
             loop = (int)(se - s) - 3;
             sen = se;
             sn = s + 2;
             while (loop >= 6)
              { if ((cn = vbp[sen[-1]][sn[1]]) <= 0) break;
                c += cn;
                stem++;
                loop -= 2;
                sen--;
                sn++; }
             if (c > e)
              { e = c;
                pos1 = s;
                pos2 = sen;
                vstem = stem; }}
          s--; }
      se--; }
  if (e > 0)
   { *varbp = (((int)(pos1-sb))<<10) + (((int)(pos2-sb))<<5) + vstem;
     return((double)(3*(vstem - 4))); }
  else
   { *varbp = 0;
     return(-12.0); }}


int *base_copy3(int *from, int *to, int n)
{ while (n-- > 0) *to++ = *from++;
  *to = TERM;
  return(to);  }


gene *nearest_trna_gene(data_set *d, gene *ts, int nt, gene *t, csw *sw)
{ int n,i,comp,mtrna,ioverlay,ilength;
  long a,b,c,e,score,thresh,psmax;
  static long proximity = 7*MINCTRNALEN/10;
  double energy;
  psmax = d->psmax;
  comp = t->comp;
  mtrna = sw->mtrna;
  ioverlay = sw->ioverlay;
  n = -1;
  energy = -INACTIVE;
  a = t->start;
  b = t->stop;
  thresh = b-a;
  if (b < a)
   { b += psmax;
     thresh += psmax;
     for (i = 0; i < nt; i++)
      { c = ts[i].start;
        e = ts[i].stop;
        if (e < c)
         { e += psmax;
           if (a > e) goto NXTW;
           if (b < c) goto NXTW;
           if (ts[i].genetype != tRNA) continue;
           if (ts[i].comp != comp) if (!mtrna) continue;
           if (ioverlay)
            { ilength = e - c;
              if ((2*thresh) > (5*ilength)) continue;
              if ((2*ilength) > (5*thresh)) continue; }
           score = (a >= c)?((b >= e)?e-a:thresh):((b >= e)?thresh:b-c);
           if (score >= proximity)
            if (ts[i].energy > energy)
              { n = i;
                energy = ts[i].energy; }
           NXTW:
           c -= psmax;
           e -= psmax; }
        if (a > e) continue;
        if (b < c) continue;
        if (ts[i].genetype != tRNA) continue;
        if (ts[i].comp != comp) if (!mtrna) continue;
        if (ioverlay)
         { ilength = e - c;
           if ((2*thresh) > (5*ilength)) continue;
           if ((2*ilength) > (5*thresh)) continue; }
        score = (a >= c)?((b >= e)?e-a:thresh):((b >= e)?thresh:b-c);
        if (score >= proximity)
            if (ts[i].energy > energy)
              { n = i;
                energy = ts[i].energy; } }
     a -= psmax;
     b -= psmax; }
  for (i = 0; i < nt; i++)
   { c = ts[i].start;
     e = ts[i].stop;
     if (e < c)
      { e += psmax;
        if (a > e) goto NXTN;
        if (b < c) goto NXTN;
        if (ts[i].genetype != tRNA) continue;
        if (ts[i].comp != comp) if (!mtrna) continue;
        if (ioverlay)
         { ilength = e - c;
           if ((2*thresh) > (5*ilength)) continue;
           if ((2*ilength) > (5*thresh)) continue; }
        score = (a >= c)?((b >= e)?e-a:thresh):((b >= e)?thresh:b-c);
        if (score >= proximity)
            if (ts[i].energy > energy)
              { n = i;
                energy = ts[i].energy; }
        NXTN:
        c -= psmax;
        e -= psmax; }
     if (a > e) continue;
     if (b < c) continue;
     if (ts[i].genetype != tRNA) continue;
     if (ts[i].comp != comp) if (!mtrna) continue;
     if (ioverlay)
      { ilength = e - c;
        if ((2*thresh) > (5*ilength)) continue;
        if ((2*ilength) > (5*thresh)) continue; }
     score = (a >= c)?((b >= e)?e-a:thresh):((b >= e)?thresh:b-c);
     if (score >= proximity)
            if (ts[i].energy > energy)
              { n = i;
                energy = ts[i].energy; } }
  if (n >= 0) return(ts + n);
  return(NULL); }


void overlap(data_set *d, gene ts[], int sort[], int n, int it, csw *sw)
{ int i,j,flag,cross,crosstoo;
  long a,b,e,f,a2,b2,e2,f2,psmax;
  char sname[100],s[100];
  flag = 0;
  cross = 0;
  psmax = d->psmax;
  a = ts[it].start;
  b = ts[it].stop;
  if (b < a)
   { a2 = a - psmax;
     b2 = b;
     b += psmax;
     cross = 1; }
  j = -1;
  while (++j < n)
   { i = sort[j];
     if (i == it) continue;
     e = ts[i].start;
     f = ts[i].stop;
     crosstoo = 0;
     if (f < e)
      { e2 = e - psmax;
        f2 = f;
        f += psmax;
        crosstoo = 1; }
     if (a <= f)
      if (b >= e)
       goto OV;
     if (crosstoo)
      if (a <= f2)
       if (b >= e2)
        goto OV;
     if (cross)
      { if (a2 <= f)
         if (b2 >= e)
          goto OV;
        if (crosstoo)
         if (a2 <= f2)
          if (b2 >= e2)
           goto OV; }
     continue;
     OV:
     if (!flag) fputc('\n',sw->f);
     name(ts+i,sname,1,sw);
     location(s,ts+i,sw,sname);
     fprintf(sw->f,"Overlap with %d: %s\n", j+1,s);
     flag = 1; }
 if (flag) fputc('\n',sw->f); }




gene *find_slot(data_set *d, gene *t, int *nts, csw *sw)
{ char s1[80],s2[80],s3[80],s4[80];
  gene *tn,*tsr;
  if (sw->comp)
   { t->stop = sw->start - t->start - 1;
     t->start = t->stop - t->nbase - t->nintron + 1;
     t->comp = 1; }
  else
   { t->start += sw->start;
     t->stop = t->start + t->nbase + t->nintron - 1;
     t->comp = 0; }
  if (!sw->linear)
   { t->start = sq(t->start);
     t->stop = sq(t->stop); }
  if (t->genetype == tRNA)
   tn = nearest_trna_gene(d,ts,*nts,t,sw);
  if (tn)
   { if (t->energy <= tn->energy) return(NULL);
     copy(tn->name,t->name);
     if (sw->verbose)
      { fprintf(stderr,"%s %s ",name(t,s1,0,sw),position(s3,t,sw));
        if (sw->energydisp) fprintf(stderr,"(%g) ",nenergy(t,sw));
        fprintf(stderr,"replacing %s %s",name(tn,s2,1,sw),
        position(s4,tn,sw));
        if (sw->energydisp) fprintf(stderr," (%g)",nenergy(tn,sw));
        fprintf(stderr,"\n"); }}
     else
      { if (*nts >= sw->space)
         { sw->space += NT;
           tsr = (gene *)realloc((void *)ts,sw->space*sizeof(gene));
           if (tsr == NULL)
            { fprintf(stderr,"No more memory to store detected genes\n");
              fprintf(stderr,"Gene lost\n");
              sw->space -= NT;
              return(NULL); }
           ts = tsr; }
        copy3cr(d->seqname,t->name,79);
        tn = ts + (*nts);
        *nts = (*nts) + 1;
        if (sw->verbose)
         { fprintf(stderr,"%s at %s",name(t,s1,0,sw),position(s2,t,sw));
           if (sw->energydisp) fprintf(stderr," (%g)",nenergy(t,sw));
           fprintf(stderr,"\n"); }}
  return(tn); }


int aatail(int *s, int *ext, csw *sw)
{ int score,e;
  static int A[6] = { 1,0,0,0,0,0 };
  static int C[6] = { 0,1,0,0,0,0 };
  if (sw->aataildiv)
   { score = 0;
     e = 0;
     if (A[s[3]])
      { score++;
        e = 3; }
     if (C[s[2]])
      { score++;
        if (!e) e = 2; }
     if (C[s[1]])
      { score++;
        if (!e) e = 1; }
     if (score < 2)
      if (A[*s]) score++;
     *ext = ++e;
     return(score); }
  else
   { score = 1;
     e = 1;
     if (C[s[1]])
      { score++;
        e = 2;
        if (C[s[2]])
         { score++;
           e = 3;
           if (A[s[3]])
            { score++;
              e = 4; }}}
     *ext = e;
     return(score); }}


int find_mt_trna(data_set *d, int *seq, int lseq,
                 int nto, csw *sw)
{ int nah,ndh,nch,nth,h,i,j,k,n,p,y,av,gcc,cgcc,catc,athresh;
  int igc,nbase,b8,b9,b48,b57,nc,na,nt,nti,nd,ndi,dposmap[32];
  int dl,tl,extastem,astem8,astem8d,ti,di,ser,tastem,tastem8,tastem8d;
  int astem,asteme,as,as8,aext,aext8,nbasefext,cloop,dloop,tloop,tc;
  int carm,cstem,darm,dstem,tarm,tstem,var,varbp,spacer1,spacer2,anticodon,ds;
  int dstemmotif,cloop7,mtxdetect;
  int *s,*sl,*s1,*s2,*s4,*sa,*sb,*sc,*se,*sf,*sg,*si;
  int *slm,*slm1,*sle,*slb,*sld,*sge;
  int *dpos,*cpos,*cend,*tpos,*tend,*apos1,*apos2,*aend1,*aend2;
  int *clooppos,*cloopend;
  unsigned int bondtype,abondtype,mabondtype,acbondtype,cbondtype;
  unsigned int agcat,cgcat,tgcat,dbondtype,dtbondtype,tbondtype;
  unsigned int r,ct[6],cm,cv,q,tendmap[63]; 
  double gcv,e,ec,ea,eas,ed,et,ev,energy,stem_energy;
  double darmthresh,tarmthresh,tthresh,dthresh,dtthresh,thresh;
  mt_trna_cloop chit[6];
  mt_trna_loop dhit[mtND+1];
  mt_trna_tloop thit[mtNTH+1];
  mt_trna_astem ahit[mtNA+1];
  gene *tn;
  static gene te =
   { "",{TERM},{TERM},NULL,0,0,0L,0L,7,7,1,2,1,4,7,5,7,0,0,0,5,0,5,7,
     tRNA,0.0,0,0,0 };
  static int cAI[6] = { 8,0,0,0,8,0 };
  static int cfCI[6] = { 0,16,0,0,16,0 };
  static int cRI[6] = { 8,0,4,0,8,0 };
  static int cTI[6] = { 0,0,0,16,16,0 };
  static int cYI[6] = { 0,8,0,4,8,0 };
  static int AI[6] = { 1,0,0,0,1,0 };
  static int CI[6] = { 0,1,0,0,1,0 };
  static int GI[6] = { 0,0,1,0,1,0 };
  static int TI[6] = { 0,0,0,1,1,0 };
  static int RI[6] = { 1,0,1,0,1,0 };
  static int YI[6] = { 0,1,0,1,1,0 };
  static int WI[6] = { 1,0,0,1,1,0 };
  static unsigned int template[6] = { 0,0,0,0,0,0 };
  static unsigned int At[6] = { 0,0,0,1,1,0 };
  static unsigned int Ct[6] = { 0,0,1,0,1,0 };
  static unsigned int Gt[6] = { 0,1,0,1,1,0 };
  static unsigned int Tt[6] = { 1,0,1,0,1,0 };
  static unsigned int cAt[6] = { 0,0,0,2,2,0 };
  static unsigned int cCt[6] = { 0,0,2,0,2,0 };
  static unsigned int cGt[6] = { 0,2,0,1,2,0 };
  static unsigned int cTt[6] = { 2,0,1,0,2,0 };
  static unsigned int aAt[6] = { 0,0,1,2,2,0 };
  static unsigned int aCt[6] = { 0,0,2,0,2,0 };
  static unsigned int aGt[6] = { 1,2,0,1,2,0 };
  static unsigned int aTt[6] = { 2,0,1,1,2,0 };
  static unsigned int dAt[6] = { 0,0,1,2,2,0 };
  static unsigned int dCt[6] = { 0,0,2,0,2,0 };
  static unsigned int dGt[6] = { 1,2,0,2,2,0 };
  static unsigned int dTt[6] = { 2,0,2,1,2,0 };
  static unsigned int clmotif[mtNCLM] =
   { 0x1321300,0x3321300,0x1323002 };
  static int dloopi[mt_DRLmaxlength+1][4] =
   { { -1 }, { -1 }, { -1 }, { -1 }, { -1 }, { -1 }, { -1 },
     { 0,2,-1 }, { 0,2,-1 }, { 0,2,3,-1 }, { 0,3,-1 }, { 0,3,-1 },
     { 0,3,4,-1 }, { 0,4,-1 }, { 0,5,-1 }, { 0,5,6,-1 }, { 0,5,6,-1 } };
  static int tloopa[12][4] =
   { { -1 }, { -1 }, { -1 }, { 0,1,-1 }, { 0,2,1,-1 }, { 4,3,2,-1 },
     { 4,3,-1 }, { 4,3,-1 }, { 4,3,-1 }, { 5,4,3,-1 }, { 5,4,-1 }, { 5,-1 } };
  static double dA[6] = { 1.0,0.0,0.0,0.0,1.0,0.0 };
  static double dT[6] = { 0.0,0.0,0.0,1.0,1.0,0.0 };
  static double C[6] = { 0.0,1.0,0.0,0.0,1.0,0.0 };
  static double G[6] = { 0.0,0.0,1.0,0.0,1.0,0.0 };
  static double T[6] = { 0.0,0.0,0.0,1.0,1.0,0.0 };
  static double AX[6] = { 0.0,-1.0,-1.0,-1.0,0.0,-1.0 };
  static double AX37[6] = { 0.0,-4.0,-1.0,-4.0,0.0,-4.0 };
  static double AXX[6] = { 0.0,-3.0,-1.5,-3.0,0.0,-3.0 };
  static double AXX37[6] = { 0.0,-4.0,-4.0,-4.0,0.0,-4.0 };
  static double AX7[6] = { 0.0,-0.7,-0.7,-0.7,0.0,-0.7 };
  static double CX[6] = { -2.0,0.0,-2.0,-1.0,0.0,-2.0 };
  static double CXX[6] = { -4.0,0.0,-4.0,-2.0,0.0,-4.0 };
  static double CX7[6] = { -0.7,0.0,-0.7,-0.7,0.0,-0.7 };
  static double TX[6] = { -1.0,-1.0,-1.0,0.0,0.0,-1.0 };
  static double TXX[6] = { -2.0,-2.0,-2.0,0.0,0.0,-2.0 };
  static double YX[6] = { -1.0,0.0,-1.0,0.0,0.0,-1.0 };
  static double tC[6] = { 0.0,0.01,0.0,0.0,0.01,0.0 };
  static double tG[6] = { 0.0,0.0,0.01,0.0,0.01,0.0 };
  static double tT[6] = { 0.0,0.0,0.0,0.01,0.01,0.0 };
  static double cA[6] = { 0.8,0.0,0.0,0.0,0.8,0.0 };
  static double cfC[6] = { 0.0,2.6,0.0,0.0,2.6,0.0 };
  static double cR[6] = { 0.8,-2.0,0.8,-0.8,0.8,-0.8 };
  static double cT[6] = { -0.8,0.0,-0.8,2.6,2.6,-0.8 };
  static double cY[6] = { -0.8,0.8,-0.8,0.8,0.8,-0.8 };
  static double loop_stab[41] =
  { 10.0,2.0,1.0,0.4,0.3,0.2,0.1,0.0,0.1,0.2,0.3,0.4,0.5,1.6,1.7,1.8,
    1.9,2.0,2.1,2.2,2.3,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,
    5.0,5.1,5.2,5.3,5.4,5.5,5.6,5.7 };
  static double bem[6][6] =
   { {  mtNOBOND, mtNOBOND, mtGABOND, mtATBOND, mtATBOND, mtNOBOND },
     {  mtNOBOND, mtNOBOND, mtGCBOND, mtNOBOND, mtGCBOND, mtNOBOND },
     {  mtGABOND, mtGCBOND, mtGGBOND, mtGTBOND, mtGCBOND, mtNOBOND },
     {  mtATBOND, mtNOBOND, mtGTBOND, mtTTBOND, mtATBOND, mtNOBOND },
     {  mtATBOND, mtGCBOND, mtGCBOND, mtATBOND, mtGCBOND, mtNOBOND },
     {  mtNOBOND, mtNOBOND, mtNOBOND, mtNOBOND, mtNOBOND, mtNOBOND } };
  static double hbem[5][5] =
   { {  0.0,0.0,0.0,mtBONDSTAB+0.5*mtATBOND,mtBONDSTAB+0.5*mtATBOND },
     {  0.0,0.0,mtBONDSTAB+0.5*mtGCBOND,0.0,mtBONDSTAB+0.5*mtGCBOND },
     {  0.0,mtBONDSTAB+0.5*mtGCBOND,0.0,mtBONDSTAB+0.5*mtGTBOND,
                                            mtBONDSTAB+0.5*mtGCBOND },
     {  mtBONDSTAB+0.5*mtATBOND,0.0,mtBONDSTAB+0.5*mtGTBOND,0.0,
                                            mtBONDSTAB+0.5*mtATBOND },
     {  mtBONDSTAB+0.5*mtATBOND,mtBONDSTAB+0.5*mtGCBOND,
                                            mtBONDSTAB+0.5*mtGCBOND,
					    mtBONDSTAB+0.5*mtATBOND,
					    mtBONDSTAB+0.5*mtGCBOND } };
  tarmthresh = sw->mttarmthresh;
  tthresh = sw->mttthresh;
  dthresh = sw->mtdthresh;
  dtthresh = sw->mtdtthresh;
  ds = sw->discrim;
  extastem = sw->extastem;
  cloop7 = sw->cloop7;
  mtxdetect = sw->mtxdetect;

  /* find cstems */

  sc = seq + sw->loffset;
  sl = seq + lseq - sw->roffset;
  h = sc[16];
  p = sc[15];
  j = sc[14];
  k = sc[13];
  n = sc[12];
  y = sc[11];
  ct[0] = cAt[h]|(cAt[p]<<4)|(cAt[j]<<8)|(cAt[k]<<12)|
	  (cAt[n]<<16)|(cAt[y]<<20);
  ct[1] = cCt[h]|(cCt[p]<<4)|(cCt[j]<<8)|(cCt[k]<<12)|
	  (cCt[n]<<16)|(cCt[y]<<20);
  ct[2] = cGt[h]|(cGt[p]<<4)|(cGt[j]<<8)|(cGt[k]<<12)|
	  (cGt[n]<<16)|(cGt[y]<<20);
  ct[3] = cTt[h]|(cTt[p]<<4)|(cTt[j]<<8)|(cTt[k]<<12)|
	  (cTt[n]<<16)|(cTt[y]<<20);
  ct[4] = 0;
  ct[5] = 0;
  for (; sc < sl; sc++)
   { p = sc[17];
     ct[0] = (ct[0] << 4) | cAt[p];
     ct[1] = (ct[1] << 4) | cCt[p];
     ct[2] = (ct[2] << 4) | cGt[p];
     ct[3] = (ct[3] << 4) | cTt[p];
     cm = (ct[sc[4]] >> 16) + (ct[sc[3]] >> 12) + (ct[sc[2]] >> 8) +
	     (ct[sc[1]] >> 4) + ct[*sc];

   /* 7 base cloop */

     cv = (cm & 0xf0);
     athresh = 12;
     nch = 0;

  /* exclude the following cloops */ 
  /* RRnnnNN, NRnnnYN */
  /* NRnnnNN with cstem < 3 Watson-Crick basepairs or equivalent */
  /* RYnnnYN */
  /* NYnnnNN with cstem < 1 Watson-Crick basepair or equivalent */
  /* NYnnnNN with cstem < 2 Watson-Crick basepairs or equivalent */
  /* unless cloop = CTnnnAN */

     if (RI[sc[6]])
      { if (RI[sc[5]]) goto CLOOP6;
        if (YI[sc[10]]) goto CLOOP6;
        if (cv < 0x60) goto CLOOP6; }
     else
      { if (YI[sc[10]])
         if (RI[sc[5]]) goto CLOOP6;
        if (cv < 0x40)
         { if (cv < 0x20) goto CLOOP6;
           if (sc[5] != Cytosine) goto CLOOP6;
           if (sc[6] != Thymine) goto CLOOP6;
           if (sc[10] != Adenine) goto CLOOP6;
            athresh = 11; }
        else if (cv < 0x70)
         { athresh = 11;
           k = cYI[sc[5]] + cTI[sc[6]] + cRI[sc[10]] + cAI[sc[11]];
           if (sc[6] == Cytosine)
            if (sc[5] == Cytosine)
             k += 16;
            else
             if (sc[5] == Thymine)
              if (sc[11] == Adenine)
               k += 16;
           if (cv == 0x40)
            { if (k < 40) goto CLOOP6;  }
           else 
            if (cv == 0x50)
             { if (k < 28) goto CLOOP6; }
            else
             { if (k < 20) goto CLOOP6;
               athresh = 9; }}
        else
         athresh = (cv < 10)?9:8; }
     chit[0].pos = sc;
     chit[0].stem = 5;
     chit[0].loop = 7;
     chit[0].looppos = sc + 5;
     chit[0].arm = 17;
     chit[0].end = sc + 17;
     chit[0].anticodon = (sc[7] << 4) + (sc[8] << 2) + sc[9];	
     if (bp[sc[-1]][sc[17]])
      { chit[1].pos = sc-1;
        chit[1].stem = 6;
        chit[1].loop = 7;
        chit[1].looppos = sc + 5;
        chit[1].arm = 19;
        chit[1].end = sc + 18;
        chit[1].anticodon = chit[0].anticodon;
        nch = 2; }
     else nch = 1;

   /* 6 base cloop */
   /* exclude cstem < 4 Watson-Crick basepairs or equivalent */
   /* exclude cloop = RRnnNN */
   /* exclude cloop = NNnnYY */

     CLOOP6:
     if (cloop7) goto CLOOPE;
     if ((cm & 0xf00) >= 0x800) 
      { if (!YI[sc[6]])
	 if (!YI[sc[5]])
          goto CLOOP8;
	if (!RI[sc[9]])
         if (!RI[sc[10]])
	  goto CLOOP8;
        se = sc + 20;
	sg = sc;
	while (sg < se)
	 { sf = sg + 5;
	   while (sf < (sg + 11)) 
	    { if (*sf == *sg)
               if (sf[1] == sg[1])
                if (sf[2] == sg[2])
                 if (sf[3] == sg[3])
                  if (sf[4] == sg[4])
	           { sb = sg + 5;
		     s = sf + 5;
		     i = 0;
		     while (sb < sf) 
		      if (*sb++ != *s++) 
		       if (++i > 1) goto NXSEG6;
                     goto CLOOPE; }
	      NXSEG6:
	      sf++; }
	   sg++; }
	chit[nch].pos = sc;
        chit[nch].stem = 5;
        chit[nch].loop = 6;
        chit[nch].looppos = sc + 5;
        chit[nch].arm = 16;
        chit[nch].end = sc + 16;
        chit[nch++].anticodon = 0;
        if (athresh > 10) athresh = 10;
        if (bp[sc[-1]][sc[16]])
         { chit[nch].pos = sc-1;
           chit[nch].stem = 6;
           chit[nch].loop = 6;
           chit[nch].looppos = sc + 5;
           chit[nch].arm = 18;
           chit[nch].end = sc + 17;
           chit[nch++].anticodon = 0; }}

   /* 8 base cloop */
   /* exclude cstem < 4 Watson-Crick basepairs or equivalent */
   /* exclude cloop = RRnnnnNN */
   /* exclude cloop = NNnnnnYY */

     CLOOP8:
     if ((cm & 0xf) >= 0x8)
      { if (!YI[sc[5]])
	 if (!YI[sc[6]])
          goto CLOOPE;
	if (!RI[sc[12]])
         if (!RI[sc[11]])
	  goto CLOOPE;
        se = sc + 20;
	sg = sc;
	while (sg < se)
	 { sf = sg + 5;
	   while (sf < (sg + 11)) 
	    { if (*sf == *sg)
               if (sf[1] == sg[1])
                if (sf[2] == sg[2])
                 if (sf[3] == sg[3])
                  if (sf[4] == sg[4])
	           { sb = sg + 5;
		     s = sf + 5;
		     i = 0;
		     while (sb < sf) 
		      if (*sb++ != *s++) 
		       if (++i > 1) goto NXSEG8;
                     goto CLOOPE; }
	      NXSEG8:
	      sf++; }
	   sg++; }
        chit[nch].pos = sc;
        chit[nch].stem = 5;
        chit[nch].loop = 8;
        chit[nch].looppos = sc + 5;
        chit[nch].arm = 18;
        chit[nch].end = sc + 18;
        chit[nch++].anticodon = 0;
        if (athresh > 10) athresh = 10;
        if (bp[sc[-1]][sc[18]])
         { chit[nch].pos = sc-1;
           chit[nch].stem = 6;
           chit[nch].loop = 8;
           chit[nch].looppos = sc + 5;
           chit[nch].arm = 20;
           chit[nch].end = sc + 19;
           chit[nch++].anticodon = 0; }}

  /* calculate carm energy */

     CLOOPE:
     if (nch < 1) continue; 
     for (nc = 0; nc < nch; nc++)
           { s1 = chit[nc].pos;
             cstem = chit[nc].stem;
             cloop = chit[nc].loop;
             s4 = s1 + cstem;
             s2 = s4 + cloop;
	     energy = (cloop == 7)?0.0:-4.0;
             energy += cY[*s4] + cT[s4[1]] + cR[s2[-2]] + cA[s2[-1]];
             if (s4[1] == Cytosine)
              if (*s4 == Cytosine)
               energy += 2.6;
              else
               if (*s4 == Thymine)
                if (s2[-1] == Adenine)
                 energy += 2.6;
             s2 += cstem;
           stem_energy = bem[*s1][*--s2];
           k = neighbour_map[*s1][*s2];
           stem_energy += neighbour_em[k][s1[1]][s2[-1]];
           bondtype = btmap[*s1][*s2];
           if (bp[*s1][*s2])
            { if (assymst[s2[1]][s1[-1]]) stem_energy += mtTERMSTAB;
              else stem_energy += send_em[*s2][*s1]; }
           else
            { if (assymst[*s2][*s1]) stem_energy += mtTERMSTAB;
              else stem_energy += send_em[s2[-1]][s1[1]]; }
           while (++s1 < s4)
          { if (!wcbp[*s1][*--s2])
             { if (!wcbp[s1[-1]][s2[1]])
    { for (j = 0; j < mtNTM; j++)
       if (*s1 == tandemid[j][1])
        if (*s2 == tandemid[j][3])
                       if (s1[-1] == tandemid[j][0])
                        if (s2[1] == tandemid[j][2])
           { stem_energy += tandem_em[j];
             break; }
      if (s1 < (s4-1))
       if (!bp[s1[1]][s2[-1]]) stem_energy -= mt3MMSTAB; }
               k = neighbour_map[*s1][*s2];
               stem_energy += (neighbour_em[k][s1[-1]][s2[1]] +
                     neighbour_em[k][s1[1]][s2[-1]]); }
           bondtype += btmap[*s1][*s2];
             stem_energy += bem[*s1][*s2]; }
      if (!bp[*--s1][*s2])
       { s1--;
         s2++; }
      if (assymst[s1[1]][s2[-1]]) stem_energy += mtTERMSTAB;
             else stem_energy += send_em[*s1][*s2];
             cgcc = bondtype & 0xf;
             if (cgcc <= 0)
              { catc = (bondtype & 0xf0) >> 4;
                if (catc < cstem) energy -= mtGCPENALTY; }
             if (cstem == 6) energy += 1.0;
      chit[nc].bondtype = bondtype;
      chit[nc].stem_energy = stem_energy;
      chit[nc].energy = energy + stem_energy; }

  /* find tarms */

  nth = 0;
  slm = sc + 61; 
  sle = sc + 57;
  sb = sc + 21;
  sg = sc + 16;
  sge = sg + 30;
  slb = sg + 32;
  template[0] = At[*slm];
  template[1] = Ct[*slm];
  template[2] = Gt[*slm];
  template[3] = Tt[*slm];
  while (--slm > sle)
   { template[0] = (template[0] << 4) | At[*slm];
     template[1] = (template[1] << 4) | Ct[*slm];
     template[2] = (template[2] << 4) | Gt[*slm];
     template[3] = (template[3] << 4) | Tt[*slm]; }
  while (slm >= sb)
   { template[0] = ((template[0] << 4) | At[*slm]) & 0xfffff;
     template[1] = ((template[1] << 4) | Ct[*slm]) & 0xfffff;
     template[2] = ((template[2] << 4) | Gt[*slm]) & 0xfffff;
     template[3] = ((template[3] << 4) | Tt[*slm]) & 0xfffff;
     sf = slm + 3;
     if (sf > sge) sf = sge;
     apos2 = slm + 5;
     si = sg;
     s = si + 4;
     r = template[*si];
     while (++si < s) r = (r >> 4) + template[*si];
     while (si <= sf)
      { if (si < slm)
          r = (r >> 4) + template[*si++];
        else
         { si++;
           r = r >> 4; }
     q = r & 0xf;
     if (slm > slb)
      { if (q < 5) continue;
        tloop = (int)(slm - si); }
     else
      { if (q < 2) continue;
        if (q < 3)
         { if (!wcbp[si[-5]][apos2[-1]]) continue;
           if (!wcbp[si[-4]][apos2[-2]]) continue;
           tloop = (int)(slm - si);
	   if (tloop > 5) continue; }
        else 
         { tloop = (int)(slm - si);
           if (q < 4)
            if (!bp[si[-4]][apos2[-2]])
             if (!bp[si[-2]][apos2[-4]])
              { if (tloop < 4) continue;
                if (si[-1] != Guanine) continue;
                if (*si != Thymine) continue;
                if (si[1] != Thymine) continue; }}}
     if (tloop < 7)
      { if (tloop < 2)
         if (tloop <= 0)
          { if (tloop <= -2)
             { if (!wcbp[si[-5]][apos2[-1]]) continue;
               if (!wcbp[si[-4]][apos2[-2]]) continue;
               tstem = 2;
               tloop += 6; }
            else
             if (bp[si[-3]][apos2[-3]])
              { tstem = 3;
                tloop += 4; }
             else
              { if (!wcbp[si[-5]][apos2[-1]]) continue;
                if (!wcbp[si[-4]][apos2[-2]]) continue;
                tstem = 2;
                tloop += 6; }}
         else
          { if (bp[si[-2]][apos2[-4]])
             { tstem = 4;
               tloop += 2; }
            else
             if (bp[si[-3]][apos2[-3]])
             { tstem = 3;
               tloop += 4; }
             else
              { if (!wcbp[si[-5]][apos2[-1]]) continue;
                if (!wcbp[si[-4]][apos2[-2]]) continue;
                tstem = 2;
                tloop += 6; }}
        else
         { if (bp[si[-1]][apos2[-5]])
            { if (q != 4) tstem = 5;
              else
               { if (bp[si[-2]][apos2[-4]]) tstem = 5;
                 else
                  { k = GI[si[-3]] + TI[si[-2]] + TI[si[-1]] + CI[*si];
                    if (k >= 2)
                     { tstem = 3;
                       tloop += 4; }
                    else tstem = 5; }}}
           else
            { if (bp[si[-2]][apos2[-4]])
               { tstem = 4;
                 tloop += 2; }
              else
               if (bp[si[-3]][apos2[-3]])
                { tstem = 3;
                  tloop += 4; }
               else
                { if (!wcbp[si[-5]][apos2[-1]]) continue;
                  if (!wcbp[si[-4]][apos2[-2]]) continue;
                  tstem = 2;
                  tloop += 6; }
	    }}
        if (tloop < 3)
         if (tstem > 3)
          { tstem--;
            tloop += 2; }}
     else
      { if (!bp[si[-1]][apos2[-5]])
         if (!bp[si[-2]][apos2[-4]])
          { tstem = 3;
            tloop += 4; }
         else
          { tstem = 4;
            tloop += 2; }
        else tstem = 5; }
     if (tloop > 17)
      if (tstem < 5)
       continue;

  /* calculate tarm energy */

     s1 = si - 5;
     tpos = s1;
     s4 = s1 + tstem;
     s2 = apos2;
     if (tt[*s1][*--s2])
      { energy = mtTSTTSTAB;
        if (tt[*++s1][*--s2])
	 { energy += mtTSTTSTAB;
	   bondtype = btmap[*s1++][*s2--]; }
	else bondtype = 0; }
     else
      { energy = 0.0;
	bondtype = 0; }

  /* calculate tstem energy */

     stem_energy = bem[*s1][*s2];
     k = neighbour_map[*s1][*s2];
     stem_energy += neighbour_em[k][s1[1]][s2[-1]];
     bondtype += btmap[*s1][*s2];
     while (++s1 < s4)
      { if (!wcbp[*s1][*--s2])
         { if (!wcbp[s1[-1]][s2[1]])
            { for (j = 0; j < mtNTM; j++)
              if (*s1 == tandemid[j][1])
               if (*s2 == tandemid[j][3])
                if (s1[-1] == tandemid[j][0])
                 if (s2[1] == tandemid[j][2])
                  { stem_energy += tandem_em[j];
                    break; }
              if (s1 < (s4-1))
               if (!bp[s1[1]][s2[-1]]) stem_energy -= mt3MMSTAB; }
           k = neighbour_map[*s1][*s2];
           stem_energy += (neighbour_em[k][s1[-1]][s2[1]] +
                     neighbour_em[k][s1[1]][s2[-1]]); }
        bondtype += btmap[*s1][*s2];
        stem_energy += bem[*s1][*s2]; }
      s1--;
      if (tloop < 4) stem_energy += ssend_em[*s1][*s2];
      else
       if (assymst[s1[1]][s2[-1]]) stem_energy += mtTERMSTAB;
       else stem_energy += send_em[*s1][*s2];

  /* compile possible tarms */

      energy += (stem_energy - mtBONDSTAB*(double)(5-tstem));
      if (energy >= tarmthresh)
       { thit[nth].pos = tpos;
         s1 = tpos + tstem;
	 s2 = apos2 - tstem;
         thit[nth].energy = energy - loop_stab[tloop] +
                  tG[s1[-1]] + tT[*s1] + tT[s1[1]] + tC[s1[2]];
         thit[nth].stem_energy = stem_energy;
         thit[nth].bondtype = bondtype;
         thit[nth].stem = tstem;
         thit[nth].loop = tloop;
	 thit[nth].end = tpos + 2*tstem + tloop;
         if (++nth >= mtNTH)
          { fprintf(stderr,"Too many mt-tstem hits\n");
            break; }
         if (tstem > 2)
	  if (tloop < 10)
           if (gt[s1[-1]][*s2])
           { thit[nth].pos = tpos;
             thit[nth].energy = energy - mtBONDSTAB - mtGTBOND -
		                loop_stab[tloop+2] +
                                tG[s1[-2]] + tT[s1[-1]] + tT[*s1] + tC[s1[1]];
             thit[nth].stem_energy = stem_energy - mtGTBOND;
             thit[nth].bondtype = bondtype - 0x100;
             thit[nth].stem = tstem - 1;
             thit[nth].loop = tloop + 2;
	     thit[nth].end = thit[nth-1].end;
             if (++nth >= mtNTH)
              { fprintf(stderr,"Too many mt-tstem hits\n");
                break; }
             if (tstem > 3)
	      if (tloop < 8)
	       if (gt[s1[-2]][s2[1]])
               { thit[nth].pos = tpos;
                 thit[nth].energy = energy - 2.0*mtBONDSTAB - 2.0*mtGTBOND -
		                loop_stab[tloop+4] + tG[s1[-3]] +
				tT[s1[-2]] + tT[s1[-1]] + tC[*s1];
                 thit[nth].stem_energy = stem_energy - 2.0*mtGTBOND;
                 thit[nth].bondtype = bondtype - 0x200;
                 thit[nth].stem = tstem - 2;
                 thit[nth].loop = tloop + 4;
	         thit[nth].end = thit[nth-1].end;
                 if (++nth >= mtNTH)
                  { fprintf(stderr,"Too many mt-tstem hits\n");
                    break; }}}
         if (tstem < 5)
          { if (tloop < 11) continue;
            if (tloop > 16) continue;
            if (!wcbp[s1[1]][s2[-2]]) continue;
            bondtype += btmap[*s1][s2[-1]] + btmap[s1[1]][s2[-2]];
            tstem += 2;
            tloop -= 4; }
         else
          { if (tloop < 9) continue;
            if (wcbp[*s1][s2[-1]])
             { if (tloop > 14) continue;
               tstem++;
               tloop -= 2;
               bondtype += btmap[*s1][s2[-1]]; }
            else
             { if (tloop < 11) continue;
               if (tloop > 16) continue;
               if (!wcbp[s1[1]][s2[-2]]) continue;
               bondtype += btmap[*s1][s2[-1]] + btmap[s1[1]][s2[-2]];
               tstem += 2;
               tloop -= 4; }}
         thit[nth].pos = tpos;
         s1 = tpos + tstem;
         thit[nth].energy = energy - loop_stab[tloop] +
                  tG[s1[-1]] + tT[*s1] + tT[s1[1]] + tC[s1[2]];
         thit[nth].stem_energy = stem_energy;
         thit[nth].bondtype = bondtype;
         thit[nth].stem = tstem;
         thit[nth].loop = tloop;
         thit[nth].end = thit[nth-1].end;
         if (++nth >= mtNTH)
          { fprintf(stderr,"Too many mt-tstem hits\n");
            break; }
         if (tloop < 9) continue;
         if (!wcbp[*s1][apos2[-tstem-1]]) continue;
         if (++tstem > 7) continue;
	 if (tloop > 14) continue;
         tloop -= 2;
         thit[nth].pos = tpos;
         s1 = tpos + tstem;
         thit[nth].energy = energy - loop_stab[tloop] +
                  tG[s1[-1]] + tT[*s1] + tT[s1[1]] + tC[s1[2]];
         thit[nth].stem_energy = stem_energy;
         thit[nth].bondtype = bondtype;
         thit[nth].stem = tstem;
         thit[nth].loop = tloop;
         thit[nth].end = thit[nth-1].end;
         if (++nth >= mtNTH)
          { fprintf(stderr,"Too many mt-tstem hits\n");
            break; }}}
     slm--; }

  /* find darms */

  ndh = 0;
  sle = sc - 4;
  slb = sc - 8;
  slm = sc - 1;
  template[0] = dAt[*slm];
  template[1] = dCt[*slm];
  template[2] = dGt[*slm];
  template[3] = dTt[*slm];
  while (--slm > sle)
   { template[0] = (template[0] << 4) | dAt[*slm];
     template[1] = (template[1] << 4) | dCt[*slm];
     template[2] = (template[2] << 4) | dGt[*slm];
     template[3] = (template[3] << 4) | dTt[*slm]; }
  slm1 = slm;
  while (slm > slb)
   { template[0] = ((template[0] << 4) | dAt[*slm]) & 0xffff;
     template[1] = ((template[1] << 4) | dCt[*slm]) & 0xffff;
     template[2] = ((template[2] << 4) | dGt[*slm]) & 0xffff;
     template[3] = ((template[3] << 4) | dTt[*slm]) & 0xffff;
     slm--;
     si = slm - 18;
     s = si + 3;
     r = template[*si];
     while (++si < s) r = (r >> 4) + template[*si];
     while (si <= slm1)
      { if (si < slm) r = (r >> 4) + template[*si++];
        else
         { r = r >> 4;
	   si++; }
        if ((q = (r & 0xf)) < 6)
         { q += (unsigned int)(TI[si[-6]] + RI[si[-5]]);
           if (q < 6) continue; }

  /* calculate darm energy */

        s1 = si - 4;
        dhit[ndh].pos = s1;
        energy = dT[s1[-2]] + dA[s1[-1]];
        dloop = (int)(slm1 - si);
        if (dloop > 2)
         if (bp[si[-1]][*slm1])
          { dstem = 4;
	    goto EC; }
        if (dloop > 0)
         if ((ggstembp[si[-2]][slm[2]]) || (gabp[si[-1]][*slm1]))
            { dstem = 3;
              dloop += 2;
              energy += mtNOBOND;
	      goto EC; }
	if (!wcbp[si[-3]][slm[3]]) continue;
	if (!gc[si[-4]][slm[4]]) continue;
        dstem = 2;
        dloop += 4;
	if (dloop > 5) energy += mtNOBOND;
        energy += mtNOBOND;

        EC:
        s2 = slm + 4;
        s4 = s1 + dstem;
        if (!wcbp[s1[1]][s2[-1]])
         if (stemterm[s1[1]][s2[-1]]) energy -= 1.0;
         else 
          if (bp[s1[1]][s2[-1]]) energy -= 1.5;
          else energy -= 2.0;

  /* calculate dstem energy */

        stem_energy = bem[*s1][*s2];
        k = neighbour_map[*s1][*s2];
        stem_energy += neighbour_em[k][s1[1]][s2[-1]];
        bondtype = btmap[*s1][*s2];
        if (bp[*s1][*s2])
         { if (assymst[s2[1]][s1[-1]]) stem_energy += mtTERMSTAB;
           else stem_energy += send_em[*s2][*s1];
           s1++;
           s2--; }
        else
         { s1++;
           s2--;
           if (assymst[s2[1]][s1[-1]]) stem_energy += mtTERMSTAB;
           else stem_energy += send_em[*s2][*s1]; }
        stem_energy += bem[*s1][*s2];
        k = neighbour_map[*s1][*s2];
        stem_energy += (neighbour_em[k][s1[-1]][s2[1]] +
                       neighbour_em[k][s1[1]][s2[-1]]);
        bondtype += btmap[*s1][*s2];
        while (++s1 < s4)
         { if (!wcbp[*s1][*--s2])
            { if (!wcbp[s1[-1]][s2[1]])
               { for (j = 0; j < mtNTM; j++)
                  if (*s1 == tandemid[j][1])
                   if (*s2 == tandemid[j][3])
                    if (s1[-1] == tandemid[j][0])
                     if (s2[1] == tandemid[j][2])
                      { stem_energy += tandem_em[j];
                        break; }
                 if (s1 < (s4-1))
                  if (!bp[s1[1]][s2[-1]]) stem_energy -= mt3MMSTAB; }
              k = neighbour_map[*s1][*s2];
              stem_energy += (neighbour_em[k][s1[-1]][s2[1]] +
                              neighbour_em[k][s1[1]][s2[-1]]); }
           bondtype += btmap[*s1][*s2];
           stem_energy += bem[*s1][*s2]; }
        if (!bp[*--s1][*s2])
         { s1--;
           s2++; }
        if (dloop < 4) stem_energy += ssend_em[*s1][*s2];
        else
         if (assymst[s1[1]][s2[-1]]) stem_energy += mtTERMSTAB;
         else stem_energy += send_em[*s1][*s2];

  /* compile possible darms */

              energy += stem_energy;
              dhit[ndh].energy = energy;
              dhit[ndh].stem_energy = stem_energy;
              dhit[ndh].bondtype = bondtype;
              dhit[ndh].stem = dstem;
              dhit[ndh].loop = dloop;
              if (++ndh >= mtND)
               { fprintf(stderr,"Too many mt-dstem hits\n");
                 break; }
	      if (dstem == 4)
	       { if (dloop >= 6)
	          if (bondtype < 0x1000)
               { s1 = si - 5;
	         s2 = slm + 5;
                 if (bp[*s1][*s2])
                  { dhit[ndh].pos = s1;
		    e = 0.5 + bem[*s1][*s2];
                    dhit[ndh].energy = energy + e;
		    if (wcbp[*s1][*s2])
		     dhit[ndh].energy += (dT[s1[-2]] + dA[s1[-1]] -
			                  dT[s1[-1]] - dA[*s1]);
                    dhit[ndh].stem_energy = stem_energy + e;
                    dhit[ndh].bondtype = bondtype + btmap[*s1][*s2];
                    dhit[ndh].stem = 5;
                    dhit[ndh].loop = dloop;
                    if (++ndh >= mtND)
                     { fprintf(stderr,"Too many mt-dstem hits\n");
                       break; }}}}
              else
               if (dloop >= 6)
               { s1 = si - 1;
	         s2 = slm1;
                 if (stemterm[*s1][*s2])
                  { dhit[ndh].pos = si - 4;
                    dhit[ndh].energy = energy;
                    dhit[ndh].stem_energy = stem_energy;
                    dhit[ndh].bondtype = bondtype;
                    dhit[ndh].stem = 4;
                    dhit[ndh].loop = dloop - 2;
                    if (++ndh >= mtND)
                     { fprintf(stderr,"Too many mt-dstem hits\n");
                       break; }}}
              if (dloop >= 4) continue;
              s1 = si - 4 + dstem - 1;
              s2 = s1 + dloop + 1;
              if (bp[*s1][*s2]) continue;
              dhit[ndh].pos = si - 4;
              dhit[ndh].energy = energy + 0.001;
              dhit[ndh].stem_energy = stem_energy;
              dhit[ndh].bondtype = bondtype;
              dhit[ndh].stem = dstem - 1;
              dhit[ndh].loop = dloop + 2;
              if (++ndh >= mtND)
               { fprintf(stderr,"Too many mt-dstem hits\n");
                 break; } }
     slm1--; }

  /* build darm exclusion map */
  /* 5' astems further from carm than */
  /* mt_DRLmaxlength must match a darm */

  for (i = 3; i <= 30; i++) dposmap[i] = 0;
  sf = sc - mt_DRLmaxlength - 1;
  sld = sf;
  if (ndh > 0)
   { s = dhit[0].pos;
     for (nd = 0; nd < ndh; nd++)
      { se = dhit[nd].pos;
        if (se < s) s = se;
	i = (int)(sc - se);
	if (dposmap[++i] < 1) dposmap[i] = 1;
	dposmap[++i] = 2;
	if (dposmap[++i] < 1) dposmap[i] = 1; }
     s -= 4;
     if (s < sf) sf = s; }
	
  /* build tarm exclusion map */
  /* 3' astems further from carm than */
  /* mt_TVRLmaxlength must match a tarm */

  for (i = 17; i <= 62; i++) tendmap[i] = 0;
  s2 = sc + mt_TVRLmaxlength + 17;
  sle = s2;
  if (nth > 0)
   { s = thit[0].end;
     for (nt = 0; nt < nth; nt++)
      { se = thit[nt].end;
	    if (se > s) s = se;
	    i = (int)(se - sc);
	    bondtype = thit[nt].bondtype;
	    if (tendmap[i])
         { if (bondtype < tendmap[i]) tendmap[i] = bondtype; }
	    else tendmap[i] = bondtype; }
     if (s > s2) s2 = s; }

  /* find astems in 3 categories: */
  /* high energy astems close to carm */
  /* high energy astems matching a high energy tarm far from carm */
  /* low energy astem matching a darm and tarm */

     nah = 0;
     sa = sc - 3;
     sg = sf - 6;
     sb = sc + 17;
     se = s2 + 6;
     template[0] = aAt[*se];
     template[1] = aCt[*se];
     template[2] = aGt[*se];
     template[3] = aTt[*se];
     while (--se > s2)
      { template[0] = (template[0] << 4) | aAt[*se];
        template[1] = (template[1] << 4) | aCt[*se];
        template[2] = (template[2] << 4) | aGt[*se];
        template[3] = (template[3] << 4) | aTt[*se]; }
     ti = (int)(se - sc);
     while (se >= sb)
      { template[0] = ((template[0] << 4) | aAt[*se]) & 0xfffffff;
        template[1] = ((template[1] << 4) | aCt[*se]) & 0xfffffff;
        template[2] = ((template[2] << 4) | aGt[*se]) & 0xfffffff;
        template[3] = ((template[3] << 4) | aTt[*se]) & 0xfffffff;
	if (tendmap[ti])
	 { nti = (tendmap[ti] < 0x2000)?1:0; }
	else
	 { if (se > sle) goto ANX;
	   nti = -1; }
        si = sg;
        r = template[*si];
        while (++si < sf) r = (r >> 4) + template[*si];
		di = (int)(sc - si);
        while (si < sa)
         { r = (r >> 4) + template[*si++];
	   if (dposmap[--di])
            { if (nti <= 0)
               { if (nti < 0)
                  if (dposmap[di] < 2) continue;
	         if ((av = (r & 0xf)) < athresh) continue; }}
	   else
	    { if (si < sld) continue;
	      if (nti < 0) continue;
              if ((av = (r & 0xf)) < athresh) continue; }
           if (nah >= mtNA)
            { fprintf(stderr,"Too many mt-astem hits\n");
              break; }

  /* predict astem length and calculate astem energy */

       s1 = si - 7;
       s2 = se + 6;
       if (bp[*s1][*s2])
        { astem = 7;
          energy = 0.0;
          ahit[nah].pos1 = s1;
          ahit[nah].pos2 = se; }
       else
        if (ggstemterm[*s1][*s2])
         { astem = 7;
           ahit[nah].pos1 = s1;
           ahit[nah].pos2 = se;
           energy = bem[*s1++][*s2--]; }
	else
         { energy = bem[*s1++][*s2--];
           if (bp[*s1][*s2])
            { astem = 6;
              ahit[nah].pos1 = s1;
              ahit[nah].pos2 = se; }
           else
            if (ggstemterm[*s1][*s2])
             { astem = 6;
               ahit[nah].pos1 = s1;
               ahit[nah].pos2 = se;
               energy += bem[*s1++][*s2--]; }
	    else
             { astem = 5;
               energy += bem[*s1++][*s2--];
               ahit[nah].pos1 = s1;
               ahit[nah].pos2 = se; }}
       ahit[nah].stem = astem;
       bondtype = btmap[*s1][*s2];
       energy += bem[*s1][*s2];
       k = neighbour_map[*s1][*s2];
       energy += neighbour_em[k][s1[1]][s2[-1]];
       energy += bem[*++s1][*--s2];
       k = neighbour_map[*s1][*s2];
       energy += (neighbour_em[k][s1[-1]][s2[1]] +
                 neighbour_em[k][s1[1]][s2[-1]]);
       bondtype += btmap[*s1][*s2];
       while (++s1 < si)
        { if (!wcbp[*s1][*--s2])
           { if (!wcbp[s1[-1]][s2[1]])
              { for (j = 0; j < mtNTM; j++)
                 if (*s1 == tandemid[j][1])
                  if (*s2 == tandemid[j][3])
                   if (s1[-1] == tandemid[j][0])
                    if (s2[1] == tandemid[j][2])
                     { energy += tandem_em[j];
                       break; }
                if (s1 < (si-1))
                 if (!bp[s1[1]][s2[-1]]) energy -= mt3MMSTAB; }
             k = neighbour_map[*s1][*s2];
             energy += (neighbour_em[k][s1[-1]][s2[1]] +
             neighbour_em[k][s1[1]][s2[-1]]); }
          bondtype += btmap[*s1][*s2];
          energy += bem[*s1][*s2]; }
       if (!bp[*--s1][*s2])
        if (!bp[*--s1][*++s2])
         if (!bp[*--s1][*++s2])
          if (!bp[*--s1][*++s2])
           goto NOST;
       if (assymst[s1[1]][s2[-1]]) energy += mtTERMSTAB;
       NOST:
       ahit[nah].energy = energy;
       ahit[nah].bondtype = bondtype;
       nah++; }
      ANX:
      se--;
      ti--; }
    if (nah <= 0) continue;

	
  /* build mttrna genes */
  /* cycle through astems first so that */
  /* GC content is only calculated once per astem */

  thresh = -INACTIVE;
  te.ps = NULL;
  for (na = 0; na < nah; na++)
      { apos2 = ahit[na].pos2;
        apos1 = ahit[na].pos1;
        astem = ahit[na].stem;
        aend1 = apos1 + astem;
        astem8 = (astem == 7)?(wcbp[apos1[-1]][apos2[7]]):0;
        asteme = 0;
        ea = ahit[na].energy;
        abondtype = ahit[na].bondtype;
	agcat = ((abondtype >> 4) + abondtype) & 0xf;

  /* GC content */

        s = apos1;
	aend2 = apos2 + astem;
	nbase = (int)(aend2 - apos1) + 1;
	igc = 0;
	while (s <= aend2)
         { k = *s++;
           if (k >= Cytosine) if (k <= Guanine) igc++; }
        gcv = 10.0*(double)igc/(double)nbase;
	if (gcv < 1.0)
         { if (gcv < 0.55) continue;
	   ea -= 0.5; }
	if (nbase > 60)
	 { if (gcv > 6.0) ea -= 2.0*(gcv - 6.0); }
	else
	 { if (gcv > 5.0) ea -= 2.0*(gcv - 5.0); }
	if (gcv > 6.6)
         { ea -= 6.0;
           if (gcv > 7.0) ea -= 6.0; }

  /* cycle through carms that fall between astem */

        nc = -1;
        while (++nc < nch)
         { cpos = chit[nc].pos;
           dloop = (int)(cpos - aend1);
           if (dloop < 3) continue;
           if (dloop > 26) continue;
           cend = chit[nc].end;
           tloop = (int)(apos2 - cend);
           if (tloop < 5) continue;
	   cloop = chit[nc].loop;
           cstem = chit[nc].stem;
           clooppos = chit[nc].looppos;
	   cloopend = clooppos + cloop;
           carm = chit[nc].arm;
           anticodon = chit[nc].anticodon;
	   cbondtype = chit[nc].bondtype;
	   acbondtype = abondtype + cbondtype;
	   cgcat = ((cbondtype >> 4) + cbondtype) & 0xf;
           ec = ea + chit[nc].energy;

  /* astem,cstem stability (GC bond count) */

         if ((abondtype & 0xf) <= 0)
          if ((cbondtype & 0xf) <= 0)
           { ec -= mtGCPENALTYD;
             if (((cbondtype & 0xf0) >> 4) >= 5) ec += 0.5; }
		

  /* anticodon to astem discriminator base match */

           astem8d = 0;
		   if (cloop == 7)
            { if (!mt_discrim[ds][anticodon][apos2[astem]])
               if (astem8)
                if (mt_discrim[ds][anticodon][apos2[8]]) astem8d = 1;
                else ec -= 3.0;
               else
                if (astem <= 6)
         { if (!mt_discrim[ds][anticodon][apos2[7]])
            if (astem == 5)
	     { if (!mt_discrim[ds][anticodon][apos2[6]]) ec -= 3.0; }
            else
             ec -= 3.0; }
        else
         ec -= 3.0; }
		
  /* build TV-replacement loop mttrna genes */

        if (tloop <= mt_TVRLmaxlength)
         { if (!sw->tvloop) goto TVN;

  /* astem termination */
  /* (only need to calculate once per astem) */

    if (!asteme)
            { asteme = 1;
              s = aend1 - 1;
              se = apos2;
              while (!bp[*s][*se])
               { if (--s <= apos1)
		  { eas = 0.0;
		    goto NOST2; }
                 se++; }
              if (!aastemterm[s[1]][se[-1]]) eas = -0.5;
      else
       { eas = 0.0;
         while (se >= apos2)
                 { s++;
                   se--;
         if (aastemterm[*s][*se]) eas += 1.0;  }}}

  /* choose darm */

        NOST2:
           energy = 94.0 + ec + eas;
           nd = -1;
           ndi = -1;
           ed = -INACTIVE;
           while (++nd < ndh)
            { dpos = dhit[nd].pos;
              spacer1 = (int)(dpos - aend1);
              if (spacer1 != 2) continue;
              dl = dhit[nd].loop;
              dstem = dhit[nd].stem;
	      if (dstem > 4) continue;
              darm = 2*dstem + dl;
              spacer2 = (int)(cpos - dpos) - darm;

  /* astem,darm,cstem interspacing */

              if (spacer2 < 1) continue;
              e = dhit[nd].energy;
              if (spacer2 > 1)
               { if (spacer2 > 2) continue;
                 if (!stembp[*cpos][cend[-1]]) continue; 
		 if (tloop > 12) e -= 2.0;
	         if ((dhit[nd].bondtype & 0xf) < 1) 
		  if ((agcat + cgcat + 1) < (cstem + astem))
                   e -= 3.0; }
              else
               if (dl > 11)
	        { if (!RI[cpos[-1]]) e -= 2.0; }
               else
                { if (cpos[-1] == Cytosine) e -= 2.0; }

  /* small,large dloop, dstem R motif */

              if (dl < 3) e -= 2.0;
	      if (dl > 12) e -= 2.0;
              if (!RI[*dpos]) e -= 1.0;

  /* darm,tloop tertiary interaction */

              k = 0;
	          di = ((dl >= 12)?3:((dl >= 9)?2:1));
              tl = (tloop >= 14)?5:((dl >= 9)?((tloop >= 10)?4:3):3);
              if (!ggstackbp[dpos[dstem+di]][cend[tl]]) 
	           { if (tl > 3)
	              { if (!ggstackbp[dpos[dstem+di]][cend[tl-1]]) e -= 1.5;
	                else k++; }
                 else 
	              if (di > 1)
	               { if (!ggstackbp[dpos[dstem+di-1]][cend[tl]]) e -= 1.5;
	                 else k++; }
		          else
                    e -= 1.5; }   
              else k++;
              if (stemterm[dpos[dstem-1]][dpos[darm-dstem]])
               { e -= 0.5;
                 if (cend[2] == dpos[dstem-2])
                  { if (bp[cend[2]][dpos[darm-dstem+1]]) k++; }
                 else
                  { if (cend[2] == dpos[darm-dstem+1]) 
                     if (bp[cend[2]][dpos[dstem-2]]) k++; }}
              else
               { if (cend[2] == dpos[dstem-1])
                { if (!bp[cend[2]][dpos[darm-dstem]]) e -= 0.5;
                  else k++; }
               else
                { if (cend[2] != dpos[darm-dstem]) e -= 0.5;
                  else
                   if (!bp[cend[2]][dpos[dstem-1]]) e -= 0.5;
                   else k++; }}
              if (cend[1] == *dpos)
               { if (!stackbp[cend[1]][dpos[darm-1]]) e -= 0.5;
                 else k++; }
              else
               { if (cend[1] != dpos[darm-1]) e -= 0.5;
                 else
                  if (!bp[cend[1]][*dpos]) e -= 0.5;
                  else k++; }

  /* darm stability */

              dstemmotif = wcbp[dpos[1]][dpos[darm-2]];
              if (spacer2 == 2)
               if ((k < 3) || (dhit[nd].bondtype > 0x200) || (!dstemmotif))
                { if (abondtype >= 0x10000) e -= 2.0;
                  if (dstem > 3) e -= 1.0; 
		          e -= 0.5; }

  /* darm tertiary interactions */

              j = 0;
              b8 = dpos[-2];
              b9 = dpos[-1];
              if (!bp[b8][dpos[dstem]]) e -= 1.0;
              else if (wcbp[b8][dpos[dstem]]) j++;
              if (!bp[b8][dpos[darm-dstem-1]]) e-= 1.0;
              else if (wcbp[b8][dpos[darm-dstem-1]]) j++;
              if (!wcbp[dpos[2]][dpos[darm-3]])
               { if (!gastembp[b8][dpos[dstem]]) e -= 2.0;
                 else if (!gastembp[b8][dpos[darm-dstem-1]]) e -= 2.0;
	         if (!ggstembp[dpos[2]][dpos[darm-3]]) e -= 1.0; }
              else j++;
              if (!bp[b9][dpos[2]])
               { if (!bp[b9][dpos[darm-3]]) e -= 1.0;
                 else j++; }
              else j++;

/* more extensive tertiary interaction between darm,tloop */

	      if (dstemmotif)
               { if (k >= 3)
                  if (bp[dpos[2]][dpos[darm-3]])
                   { if (b8 != Thymine) e += 0.5;
                     if (dl > 3)
                      if (bp[dpos[dstem+2]][cend[tl+1]]) e += 0.7;
                      else
                       if (gabp[dpos[dstem+2]][cend[tl+1]]) e += 0.5;
		     if (tloop >= 6)
		      if (spacer2 < 2)
                       if (dl >= 3)
		        { di = (dl > 11)?2:1;
                          if (bp[dpos[dstem+di]][cend[tl]])
		           { if (chit[nc].stem_energy > -4.8) e += 0.5;
                             if (wcbp[dpos[dstem+di]][cend[tl]])
                              if (gcv > 1.2)
                               if (clooppos[1] == Thymine)
			        if (cbondtype < 0x200)
			         if ((cbondtype & 0xf) > 0)
                                  if (abondtype < 0x2000)
				   { e += 1.5;
				     if (dl > 3)
                                      if (wcbp[dpos[dstem+di+1]][cend[tl+1]])
                                       e += 1.0; }}}}
                 if (j >= 4) e += 0.25; }
              if (e > ed)
               { ed = e;
                 ndi = nd;
	         ti = k; }}
      if (ndi < 0) goto TVN;
      energy += ed;
      dpos = dhit[ndi].pos;
      dstem = dhit[ndi].stem;
      dl = dhit[ndi].loop;
      darm = 2*dstem + dl;
      dbondtype = dhit[ndi].bondtype;
      spacer2 = (int)(cpos - dpos) - darm;
      spacer1 = (int)(dpos - aend1);
      b8 = *aend1;
      b9 = aend1[1];

  /* false positive suppression */

      if (dloop < 15) energy -= 2.0;
      if (cstem > 5) energy -= 1.0;
      if (tloop < 6) energy -= 1.0;
      if (tloop > 12)
       { energy -= 1.0;
         if (agcat < 6) energy -= 2.0;
         if (tloop > 15) energy -= 2.5; }
      if (!stackbp[*dpos][dpos[darm-1]]) energy -= 1.0;
      if (dstem < 4)
       if (gcv > 1.2)
        if ((dbondtype & 0xf0f) == 0) energy -= 1.5;
     if (b8 != Thymine)
      { if (dl < 4)
	     if (abondtype > 0x10000)
           energy -= 1.5;
        if (b8 == Adenine) 
         if (YI[cloopend[-2]])
          energy -= 1.0; }  
      if (dl > 10)
       { if (tloop < 7) energy -= 2.0;
         if (spacer2 > 1) energy -= 2.0;
	 if (dhit[ndi].stem_energy < -3.4) energy -= 2.0; }
      if (gcv < 2.0)
       if (dbondtype > 0x10000) energy -= 2.0;
      if ((cbondtype & 0xf) < 1)
       if (abondtype > 0x100)
        { if (cgcat < 4)
           energy -= 1.5;
	  if (!wcbp[dpos[2]][dpos[darm-3]]) energy -= 1.0; }
      if (b8 != Thymine)
       if ((clooppos[1] != Thymine) ||
           (*clooppos != Cytosine))
        if (dl > 3)
         if (dbondtype > 0x10000)
           energy -= 1.0;
      if (!RI[cend[1]])
       if (b9 != Guanine) energy -= 1.0;
       else energy -= 0.5;
      if (b9 == Guanine)
       { if (!RI[*cend]) energy -= 1.0;
         if (spacer2 != 1) energy -= 3.0;
         else
          { tl = (tloop >= 14)?5:((dl >= 9)?((tloop >= 7)?4:3):3);
            s = dpos + dstem;
            if (!wcbp[s[1]][cend[tl]])
	     { energy -= 2.5;
               if (dl >= 5)
                if (chit[nc].energy > 2.0)
                 if (wcbp[s[2]][cend[tl]])
                  if (wcbp[s[3]][cend[tl+1]])
                   energy += 6.0; }
	    else
	     if (b8 == Thymine)
             if (dl >= 5)
              if (chit[nc].energy > 2.0)
               if (wcbp[s[2]][cend[tl+1]])
                   energy += 3.5; }}
      else if (b9 != Adenine) energy -= 3.0;
      if (b8 != Thymine)
       if (b8 == Guanine)
        { if (!RI[dpos[dstem]]) energy -= 1.0;
          else
           if (RI[dpos[darm-dstem-1]]) energy += 2.0; }
       else energy -= 1.0;

   /* carm termination */

       if (assymst[cend[-1]][*cpos]) energy += 1.0;

   /* CTnnnAA cloop motif */

       energy += CX7[*clooppos] + AX7[cloopend[-2]];
       if (clooppos[1] == Cytosine) energy -= 2.0;

   /* NNnnnAA cloop motif */

       if (cloopend[-2] == Adenine)
        if (cloopend[-1] == Adenine)
             if (spacer1 == 2)
              if (dbondtype < 0x1000)
               { if (abondtype < 0x100) energy += 1.0;
                 else
		  if (cbondtype < 0x100) energy += 1.0; }

  /* global stem damage level */

           bondtype = acbondtype + dbondtype;
           i = (int)((bondtype >> 16) & 0xf);
           j = (int)((bondtype >> 12) & 0xf);
           k = (int)((bondtype >> 8) & 0xf);
	   if (k > 0)
	    if (i > 0)
             { k += (i + j);
               if (k > 5) energy -= 1.0*(double)(k - 5); }

  /* global stem stability (GC bond count) */

      gcc = bondtype & 0xf;
      if (gcc < 2) 
       { if (ti >= 2)
          { if (cbondtype < 0x100)
             if ((cbondtype & 0xf) > 0) goto NGCC1;
            if (ti >= 3)
             if (cgcat >= 4)
              { if ((cbondtype & 0xf) > 0) goto NGCC1;
                if (cbondtype < 0x100) goto NGCC2; }}
	  energy -= (double)(3 - gcc);
	 NGCC2:
         if (gcc < 1)
          { if (agcat < 5) energy -= 2.0;
	    if (bondtype > 0x10000) energy -= 1.5; }}
	 NGCC1:


  /* global stability */
  /* (stem stability,dloop-tloop tertiary interaction,dloop size) */

	   if (abondtype > 0x1000)
            if (ti < 3)
	     { if (chit[nc].stem_energy < -6.0)
                energy -= 1.5;
               if (dl > 9)
		if (((dbondtype + cbondtype) & 0xf) < 1)
                 energy -= 1.0; }

  /* tloop,dloop tertiary interaction */
  /* (alternative dloop position) */

	   if (bondtype < 0x1000)
            if (b8 == Thymine)
             if (RI[b9])
              if (dl > 4)
               if (!bp[cend[3]][dpos[dstem+1]])
                if (bp[cend[3]][dpos[dstem+2]])
                 energy += 0.5;


  /* "near perfect" TV-loop mttRNA: */
  /* darm-tloop tertiary interaction,low global stem damage, */
  /* TR motif at b8-9, good astem,darm,carm interspacing */

           if (ti >= 2)
            if (agcat >= 6)
             if (cbondtype < 0x100)
              if (dbondtype < 0x100)
               if (RI[b9])
               if (b8 == Thymine)
                if ((abondtype & 0xf) > 0)
                 if ((dbondtype & 0xf) > 0)
                  if (spacer1 == 2)
                   if (spacer2 == 1)
                    energy += 1.5;

  /* find exceptions */

           if (energy < dthresh)
            { if (!mtxdetect) goto TVN;
	      if (energy < (thresh - 7.0)) goto TVN;
	      if (energy < (dthresh - 7.0)) goto TVN; 
	      if (nbase > 68) goto TVN;
              if (abondtype > 0x20100) goto TVN;
	      if (dl > 9) 
               { if (dl > 10) goto TVN;
                 if (dstem < 4) goto TVN;
                 if (dbondtype > 0x100) goto TVN; }
	      if (dstem > 4) goto TVN;
              if (b9 != Adenine) 
               { if (b9 != Guanine) goto TVN;
		 if (cbondtype > 0x100) goto TVN;
		 if (dbondtype > 0x200) goto TVN; }
	     if (cloop != 7) goto TVN;
         if (YI[cloopend[-2]]) goto TVN;
              if (b8 == Thymine)
                { if (apos2[-1] == Thymine)
                   if (apos2[-2] == Thymine)
                    if (tloop < 8)
                     if (tt[aend1[-1]][*apos2])
                      if (wcbp[dpos[2]][dpos[darm-3]])
                       if (((dbondtype + cbondtype) & 0xf) > 0)
			energy += 3.0; }
               else
                if (b8 == Adenine)
                 { if (apos2[-1] == Adenine)
                    if (apos2[-2] == Adenine)
		     { if (assymat[aend1[-1]][*apos2])
                        if (assymat[apos2[1]][aend1[-2]]) energy += 2.0;
		       if (agcat >= 5)
		        if (cgcat >= 4)
                         if (dbondtype < 0x100)
		          if (at[aend1[-1]][*apos2])
                           if (at[apos2[1]][aend1[-2]])
                            energy += 1.0; }
                   if (ti >= 3)
                    if (cgcat >= 4)
                     if (agcat >= 4)
                      if ((cbondtype & 0xf) > 0)
                       if ((abondtype & 0xf) > 1)
                      if (dbondtype < 0x200)
                       if (wcbp[dpos[1]][dpos[darm-2]])
                        if (clooppos[1] == Thymine)
                         if (YI[*clooppos])
                          if (RI[cloopend[-2]])
                           if (RI[cloopend[-1]])
                            energy += 5.0; }
	      if (bondtype < 0x100)
               { if (spacer2 == 1)
	          if (*clooppos == Cytosine)
		    if (clooppos[1] == Thymine)
                     if (cloopend[-2] == Adenine)
                      if (cloopend[-1] == Adenine)
                       energy += 2.0; }
              else
               { if (spacer2 == 1)
                     { if (b8 == Thymine)
		        if (dl > 3)
                         if (dbondtype < 0x200)
                         { if (cbondtype < 0x100)
                            { if (!bp[dpos[dstem+1]][cend[3]])
                               if (bp[dpos[dstem+1]][cend[4]])
                                energy += 2.0;
                              if (dbondtype < 0x100)
                               if (abondtype < 0x20000)
                                if (ti >= 2)
                                 if (dstem >= 3)
                                  if (tloop < 13)
                                   if ((cbondtype & 0xf) > 0)
                                    energy += 4.0; }}
                         else
                          if (dstem > 3)
                          if (dbondtype < 0x300)
			   { if (bondtype < 0x10000)
                              if (ti >= 3)
                               if ((acbondtype & 0xf) > 0)
                                if (wcbp[dpos[2]][dpos[darm-3]])
                                 energy += 4.0; }
                       if (tloop < 8)
		        { if (dbondtype < 0x200)
			   { if (cbondtype < 0x100)
                              if (ti >= 2)
                               {
                               if (wcbp[dpos[dstem+1]][cend[3]])
			        {
                                if (b8 == Thymine)
			         if (abondtype < 0x3000)
                                  energy += 5.0;
                                if (agcat >= 5)
			         if (gcv > 1.2)
				  if (RI[cloopend[-1]])
                                   energy += 7.0;
			        }
		               if (dbondtype < 0x100)
                                if (agcat >= 6)
                                 if (YI[*clooppos])
                                  if (clooppos[1] == Thymine)
                                   if (RI[cloopend[-2]])
                                    if (RI[cloopend[-1]])
                                     energy += 2.0; 
                                }
			     if (cbondtype < 0x300)
			      if (ti >= 3)
                               if (abondtype < 0x2000)
				if ((dbondtype & 0xf) > 0)
				 if ((acbondtype & 0xf) > 0)
                                  if (ahit[na].energy >= -7.0)
                                   if (dstem >= 4)
                                     energy += 3.0; }
			  if (dbondtype < 0x300)
                           if (cgcat >= 4)
                            if (abondtype < 0x2000)
                             if (ahit[na].energy >= -7.0)
                              if (cbondtype < 0x10000)
                               if ((cbondtype & 0xf) > 0)
                                if (cstem < 6)
			         if (ti >= 3)
			          energy += 4.0; }}
                  if (tloop > 8)
                   if (agcat >= 6)
                    if (cbondtype < 0x100)
                     if ((cbondtype & 0xf) > 0)
                      if (b8 == Thymine)
                       if (wcbp[dpos[dstem+1]][cend[3]])
                        if (wcbp[dpos[1]][dpos[darm-2]])
                         energy += 7.0; }
              if (dbondtype < 0x100)
               if (cgcat >= 4)
                if (agcat >= 5)
                        if (wcbp[dpos[1]][dpos[darm-2]])
                               if ((cbondtype & 0xf) > 0)
                               if ((abondtype & 0xf) > 0)
                               if ((dbondtype & 0xf) > 0)
				       energy += 0.5;
	      if (cbondtype < 0x100)
	       if (dbondtype < 0x200)
                if (agcat >= 5)
                 if (b8 == Thymine)
                  if (tloop < 8)
                   if (wcbp[dpos[1]][dpos[darm-2]])
                    if (wcbp[dpos[2]][dpos[darm-3]])
                     if ((cbondtype & 0xf) > 0)
                      if ((abondtype & 0xf) > 0)
                       if ((dbondtype & 0xf) > 0)
                        if (clooppos[1] == Thymine)
                         if (YI[*clooppos])
                          if (RI[cloopend[-2]])
                           energy += 3.0;
	      if (energy < dthresh) goto TVN;
              energy -= (0.9*(energy - dthresh) + 5.0); }

  /* remember fully formed TV-loop replacement mttRNA gene */
  /* if threshold reached */

         if (energy < thresh) goto TVN;
         te.energy = energy;
         thresh = energy;
         te.ps = apos1;
         te.dstem = dstem;
         te.dloop = dl;
         te.spacer1 = spacer1;
         te.spacer2 = spacer2;
         te.cstem = cstem;
         te.cloop = cloop;
         k = astem + spacer1 + darm + spacer2;
         te.anticodon = k + cstem + 2;
         te.nintron = 0;
         te.intron = 0;
         te.var = 0;
         te.varbp = 0;
         te.tstem = 0;
         te.tloop = tloop;
         te.nbase =  k + carm + tloop;
	 tastem = astem;
	 tastem8 = astem8;
	 tastem8d = astem8d;
	
  /* build D-replacement loop mttrna genes */

          TVN:
          if (tloop < 10) continue; }
          if (dloop > mt_DRLmaxlength) goto DN;
	  if (gcv < 1.2) goto DN;
          energy = 91.0 + ec;

  /* CCnnnAA cloop */

          if (clooppos[1] == Cytosine)
            { if (*clooppos != Cytosine) goto DN;
              if (cloopend[-2] != Adenine) goto DN;
              if (cloopend[-1] != Adenine) goto DN;
	      energy -= 1.0; }

  /* choose tarm */

           nt = -1;
           nti = -1;
           et = -INACTIVE;
           while (++nt < nth)
            { tl = thit[nt].loop;
              if (tl > 11) continue;
	      if (thit[nt].end != apos2) continue;
              tpos = thit[nt].pos;
              tstem = thit[nt].stem;

  /* var loop (3-7 bases long) */

              var = (int)(tpos - cend);
              if (var < 3) continue;
              e = thit[nt].energy;
              if (var > 5)
               { if (var > 7) continue;
		 if (tl < 7) continue;
		 e -= 1.0;
                 if ((dloop < 10) || (tstem < 4))
                  e -= 2.0*(double)(var - 5); }

  /* tloop RA or RG motif */

              s = tpos + tstem;
              k = 0;
	      n = 0;
              i = 0;
              while ((j = tloopa[tl][i++]) >= 0)
               if (s[j] == Adenine)
                { k = 1;
                  if (dloop >= 3)
                   if (tl > 3)
                    { b57 = s[j-1];
                      if (RI[b57] || (tl < 5))
                       { if (bp[b57][aend1[0]])
                          { e += 1.5;
                            n = 1;
                            break; }
                         if (bp[b57][aend1[1]])
                          { e += 1.5;
                            n = 1;
                            break; }
                         if (dloop > 10)
                          if (bp[b57][aend1[2]])
                           { e += 1.5;
                             n = 1;
                             break; }}}}
              if (!k)
               { i = 0;
                 while ((j = tloopa[tl][i++]) >= 0)
                  if (s[j] == Guanine)
                   if (RI[s[j-1]])
                    { k = 1;
                      break; }
                 if ( j < 0) e -= ((tl > 5)?2.0:1.0); }

  /* tertiary interaction between tloop and start of dloop */

              ti = (tl > 5)?1:((dloop > 5)?1:0);
              di = (dloop > 5)?2:1;
              if (stackbp[aend1[di]][s[ti]]) e += 1.0;

  /* tloop GTTC motif */

              i = (s[-1] == Guanine)?1:0;
              if (tl >= 5)
               { ti = i + TI[*s] + TI[s[1]] + CI[s[2]];
                 if (n)
                  if (!i)
                   if (TI[*s])
                    if (TI[s[1]])
                     if (AI[s[2]])
                      if (tl >= 7)
                       ti++;
                 if ((i > 0) || (ti >= 3)) e += (double)ti; }
              else
               { ti = i + TI[*s] + TI[s[1]];
                 if ((i > 0) || (ti >= 2)) e += (double)ti; }
              if (e > et)
               { et = e;
                 nti = nt;
                 tc = k; }}

        if (nti < 0) goto DN;
           energy += et;
           tpos = thit[nti].pos;
           tstem = thit[nti].stem;
           tl = thit[nti].loop;
           tbondtype = thit[nti].bondtype;
           var = (int)(tpos - cend);

  /* tertiary interaction between b48(=tpos[-1]) and dloop */

           b48 = tpos[-1];
           if (dloop <= 7)
            { if (YI[b48]) tc++;
              else energy -= 1.0; }
           else
            { i = 0;
              while ((j = dloopi[dloop][i++]) >= 0)
              if (assymagbp[b48][aend1[j]])
               { tc++;
                 break; }
              if (j < 0) energy -= 1.0; }

   /* large dloop, large tloop */

           if (dloop > 7)
	    { if (tl >= 6)
               if (tc < 2)
                energy -= 2.0;
	      if (tstem < 3)
                energy -= 1.0; }

   /* carm termination */

           s = cpos - 1;
           se = cend;
           if (cstem > 5)
            { s++;
              se--; }
           if (!stackbp[*s][*se]) energy -= 1.0;
           se = cpos - 3;
           if (!bp[cend[-1]][*cpos])
            { if (assymst[cend[-1]][*cpos])
               { if (dloop < 5) se++;
                 energy += 1.5; }
              else if (dloop < 13) se++; }
           else
            { if (cstem > 5) { if (dloop < 13) se++; }
              else if (dloop < 5) se++; }

   /* tertiary interaction between tloop and dloop near carm */

    s = tpos + tstem;
    if (tl >= 5)
     { ti = (tl >= 10)?4:((tl >= 7)?3:2);
       b57 = s[ti];
       if (!gabp[*se][b57]) energy -= 2.0;
       else
        { k = (var > 3)?2:((var > 1)?1:0);
          if (bp[cend[k]][b57]) energy += 1.0; }}

  /* R motif at end of tstem */

    if (!RI[s[-1]]) energy -= 2.0;

  /* large tloop */

    if (tl > 9) 
     if (tbondtype > 0x200) energy -= 2.0;

  /* dloop,var,tloop T repeat motif */
  /* present in some nematode D-loop replacement tRNA-Ser genes */

    if (dloop >= 4)
     { k = 1;
       se = aend1;
       while (se < cpos) if (*se++ == Thymine) k++;
       if (k >= dloop)
        { if (var >= 3)
           { se = cend;
             while (se < tpos) if (*se++ == Thymine) k++;
             if (k >= (var + dloop))
              { energy += 3.0;
                se = s + ((tl > 5)?5:tl);
                while (s < se) if (*s++ != Thymine) break;
                if (s >= se) energy += 5.5; }}}}

  /* astem stability  */

    if (ea < -6.1)
     if (tl > 4)
      { if (*s == Thymine)
         if (s[-1] == Guanine)
          if (s[1] == Thymine)
           goto NASI;
        if (ea > -8.3)
         if (*clooppos == Cytosine)
          if (clooppos[1] == Thymine)
           if (cloopend[-2] == Adenine)
            if (cloopend[-1] == Adenine)
             goto NASI;
        energy -= 3.0; }
    NASI:

  /* cstem stability (GC bond count) */

    bondtype = acbondtype + tbondtype;
    if ((cbondtype & 0xf) < 1)
     if ((bondtype & 0xf) < 3) energy -= 1.0;

  /* cloop CTnnnAA motif */

    if (bondtype >= 0x400)
     energy += CX[*clooppos] + TX[clooppos[1]] +
               AXX[cloopend[-1]] + AXX37[cloopend[-2]];
    else
     energy += CX[*clooppos] + TX[clooppos[1]] +
               AX[cloopend[-1]] + AX37[cloopend[-2]];

  /* large dloop */

    if (dloop >= 9)
     { k = tloop - dloop - 4;
       if (k < 0)
        if (bondtype >= 0x1000) energy += (double)k;
         if (dloop >= 12)
          { if (dloop >= 14) energy -= 2.0;
            else
             if (tstem < 6) energy -= ((dloop >= 13)?2.0:1.0); }}

  /* small dloop, small tarm */

    if (dloop <= 10)
     if (tstem < 3) 
      if (ea > -2.6)
        if (tl <= 7)
          if (cgcat >= 4)
            if (gc[*tpos][apos2[-1]])
              if (gc[tpos[1]][apos2[-2]])
               if (gcv > 1.2)
                if ((abondtype & 0xf) > 0)
                  if ((cbondtype & 0xf) > 0)
		   energy += (4.5 + (mtBONDSTAB - 0.5)*(double)(5 - tstem));

  /* global stem damage level */

           i = (int)((bondtype >> 16) & 0xf);
           j = (int)((bondtype >> 12) & 0xf) + i;
           k = (int)((bondtype >> 8) & 0xf);
	   if (tstem > 3)
	    { if ((k > 0) || (tl > 9))
	       if ((j > 0) || (k > 5))
                { n = j + k;
                  if ((s[-1] != Guanine) || (*s != Thymine) || 
                    (s[1] != Thymine) || (tstem < 5))
                   if (n > 4) energy -= 2.0*(double)(n - 4); }}
	   else
	    { n = j + k;
              if (n > 3) energy -= 2.0*(double)(n - 3); }

  /* long tstem with tloop GTT motif */

   if (s[-1] == Guanine)
    if (*s == Thymine)
     if (s[1] == Thymine)
      if (tstem >= 6)
       if (tbondtype < 0x100)
        energy += 1.5;

  /* find exceptions */

           if (energy < tthresh)
            { if (!mtxdetect) goto DN;
              if (energy < (thresh - 13.5)) goto DN;
              if (energy < (tthresh - 13.5)) goto DN;
              if (k > 1)
	       { if (i > 2) goto DN;
                 if (k > 4)
	          if (i > 1) goto DN; }
	      if (nbase > 70) goto DN;
	      if (var > 4)
	       { if (var > 5) goto DN;
	         if (var > tl) goto DN; }
	      if (tstem < 4)
               if ((agcat + cgcat + 2) < (astem + cstem)) 
		goto DN;
              if (tl > 9) goto DN; 
              if (dloop > 13) goto DN;
	      if (!YI[*clooppos]) goto DN;
	      if ((abondtype & 0xf) < 2)
	       { if ((abondtype & 0xf) < 1) goto DN;
                 if (cbondtype > 0x200)
                  if (tbondtype > 0x100)
                   if (abondtype > 0x200)
                    goto DN; }
	      if ((tbondtype & 0xf) < 1)
	       { if ((acbondtype & 0xf) < 1) goto DN;
                 if (acbondtype > 0x200) goto DN; }
              if ((dloop + 19) < tloop) goto DN;
              if (gcv > 5.5) goto DN;
              tgcat = ((tbondtype >> 4) + tbondtype) & 0xf;
	      if ((tgcat + 2) < tstem) goto DN;
	      if (cloop != 7) goto DN;
              if (bp[*cpos][cend[-1]])
                  if (bp[cpos[-1]][*cend])
                   if (bp[cpos[-2]][cend[1]])
		    energy += 2.0;
              if (bondtype < 0x20000)
                  if (thit[nti].stem_energy > -4.6)
                   if (tstem >= 4)
 	            if ((tstem >= 5) || (s[-1] == Guanine))
                     if (stackbp[cpos[1]][cend[-2]])
                      if (stackbp[*cpos][cend[-1]])
                       if (stackbp[cpos[-1]][*cend])
                        { energy += 1.5;
	                  if (s[-1] == Guanine)
                           if (*s == Thymine)
                            if (s[1] == Thymine)
                             energy += 1.0; 
                          if (agcat >= 6) energy += 0.5; }
               if (tc > 0)
                 if (tstem >= 5)
                   if (var < 6)
                    if (var > 2)
		     { if (acbondtype < 0x100) energy += 5.0;
                       else
			if ((abondtype + tbondtype) < 0x100)
			 energy += 3.0;
		        else
                         if (cloopend[-2] == Thymine)
                         if (cloopend[-1] == Thymine)
                         if (dloop > 7)
			 if (tbondtype < 0x100)
                          if (!tt[*tpos][apos2[-1]])
                          if ((agcat+cgcat) >= 10)
                           energy += 13.5; }
	         if (s[-1] == Guanine)
                  if (*s == Thymine)
                   if (s[1] == Thymine)
		    if ((tstem >= 5) || (s[2] == Cytosine))
		     { energy += 1.5;
		       if (tstem >= 5)
                        if (tbondtype < 0x1000)
		         if (s[2] == Cytosine)
                          { if (abondtype < 0x10000)
                             { if (*clooppos == Cytosine)
                                if (clooppos[1] == Thymine)
                                 if (cloopend[-2] == Adenine)
                                  if (cloopend[-1] == Adenine)
                                   energy += 3.0;
			       if (tbondtype < 0x200)
			        if (bondtype < 0x10000)
			         if (tl == 7)
			          if (s[4] == Adenine)
			           energy += 4.0; }}
			 else
			  if (tbondtype < 0x200)
                           if ((tbondtype & 0xf) >= 2)
                             if (*clooppos == Cytosine)
                              if (clooppos[1] == Thymine)
                                 if (cloopend[-2] == Adenine)
                                  if (cloopend[-1] == Adenine)
                                 energy += 1.0; }
		 if (tstem >= 4)
                  if (tbondtype < 0x100)
	           if (cbondtype < 0x200)
		    if (agcat >= 5) energy += 1.5;
	         if (energy > tthresh) energy = tthresh; 
              if (ea > -1.8) energy += 3.0;
              else
	       if (abondtype < 0x60) energy += 1.5;
	       else
		if (acbondtype < 0x200) energy += 0.75;
              if (*clooppos == Cytosine)
                                 if (cloopend[-2] == Adenine)
                                  if (cloopend[-1] == Adenine)
		 { if (tstem >= 5)
                    if (tbondtype < 0x100)
                     if (clooppos[1] == Thymine)
                      { energy += 3.0;
                        if (tstem >= 6) energy += 1.0; }
                     else
                      if (clooppos[1] == Cytosine)
                       energy += 1.0;
		   if (tc >= 2)
		    if (clooppos[1] == Thymine)
                     if (bondtype < 0x1000)
                      if (tstem >= 4)
                        if (var < 6)
                         if (var > 2)
		          energy += 3.0; }
              if (cbondtype < 0x100)
               if (agcat >= 5)
                if (tc > 0)
                 if (clooppos[1] == Thymine)
                  if (YI[*clooppos])
                   if (RI[cloopend[-2]])
                    if (RI[cloopend[-1]])
                     if (tbondtype < 0x100)
                      energy += 4.0;
	             else
                      if (agcat >= 6)
                       if ((tgcat + 1) >= tstem)
                        if (tstem >= 4)
                         energy += 4.0;
              if (bondtype < 0x1000) 
               { energy += 0.5;
                 if (bondtype < 0x200) energy += 0.75; } 
	      if (energy < tthresh) goto DN;
              energy -= (3.0 + 0.9*(energy - tthresh)); }

  /* mammalian cloop motif constraint */

	   if (ds == MAMMAL_MT)
	    { s1 = clooppos;
	      s2 = s1 + cloop;
	      r = *s1++;
	      while (s1 < s2) r = (r << 4) + *s1++;
              if (r != clmotif[0])
               if (r != clmotif[1])
                if (r != clmotif[2])
                 energy -= 5.0; }

  /* remember fully formed D-loop replacement mttRNA gene */
  /* if threshold reached */

              if (energy < thresh) goto DN;
	      te.energy = energy;
              thresh = energy;
              te.ps = apos1;
              te.spacer1 = 0;
              te.spacer2 = 0;
              te.dstem = 0;
              te.dloop = dloop;
              te.cstem = cstem;
              te.cloop = cloop;
              te.anticodon = astem + dloop + cstem + 2;
              te.nintron = 0;
              te.intron = 0;
              te.var = var;
              te.varbp = 0;
              te.tstem = tstem;
              te.tloop = tl;
              te.nbase = astem + dloop + carm + var +
                  2*tstem + tl;
	      tastem = astem;
	      tastem8 = astem8;
	      tastem8d = astem8d;

  /* build fully formed mttRNA genes */

       DN:
       if (dloop < 10) continue;

  /* choose tarm */

          nt = -1;
          nti = -1;
          et = -INACTIVE;
          while (++nt < nth)
	   { tend = thit[nt].end;
             if (tend != apos2) continue;
             e = thit[nt].energy;
	     tpos = thit[nt].pos;
             tstem = thit[nt].stem;

  /* GT motif on tloop */

	     s = tpos + tstem;
	     if (*s == Thymine)
	      if (s[-1] == Guanine)
	       if (tstem >= 5)
	        if (!stackbp[*tpos][tend[-1]])
                 { e += 0.5;
	           if (!bp[tpos[1]][tend[-2]]) e += 0.5; }

  /* large var loop */

             var = (int)(tpos - cend);
             if (var > 5)
              { ev = (double)(var - 5);
		if (tstem < 5) e -= 3.0*ev;
		else e -= (0.5 + 0.5*ev);

  /* allow large var loop if tarm looks nuclear */
  /* (GTTC motif, very large var loop base-pairing) */

                if (var > 9)
                 { if ((thit[nt].bondtype & 0xf) < 1) e -= 1.0;
                   e -= (0.25*(double)(var - 8));
		   if (*s == Thymine)
		    if (s[-1] == Guanine)
		     if (s[1] == Thymine)
		      if (s[2] == Cytosine)
		       e += 4.0;
                   if (var > 17)
                    { if (var > 25) continue;
                      e += 0.5*vloop_stability(cend,var,&varbp); }}}

  /* small var loop */

             if (var < 3)
              { if (tstem > 5)
	         if (s[-1] != Guanine)
		  e -= 0.5;
                if (var < 2)
	         { if (var < 1)
                    { if (var < 0) continue;
                      if (tstem < 4)
                       if (thit[nt].stem_energy < -4.0)
                        continue; }
                   e -= 3.0; }}
             if (e > et)
              { et = e;
                nti = nt; }}

          if (nti < 0) continue;
          tpos = thit[nti].pos;
          tstem = thit[nti].stem;
          tl = thit[nti].loop;
          tarm = 2*tstem + tl;
          var = (int)(tpos - cend);
          b48 = tpos[-1];
          tbondtype = thit[nti].bondtype;
	  bondtype = acbondtype + tbondtype;
	  ti = (int)(((bondtype >> 16) & 0xf) + ((bondtype >> 12) & 0xf) +
               ((bondtype >> 8) & 0xf));

  /* choose darm */

          nd = -1;
          ndi = -1;
          ed = -INACTIVE;
          while (++nd < ndh)
           { dl = dhit[nd].loop;
             dstem = dhit[nd].stem;
             darm = 2*dstem + dl;
             dpos = dhit[nd].pos;
             e = dhit[nd].energy;

  /* spacing between astem,darm,carm */

             spacer1 = (int)(dpos - aend1);
             spacer2 = (int)(cpos - dpos) - darm;
             if (spacer1 < 2)
              { if (spacer1 < 1) continue;
		if (dstem < 3) continue;
                if (dl > 12) e -= 2.0;
		if (astem < 7) e -= 1.0;
                if (spacer2 != 2)
                  { if (spacer2 < 1) continue;
                    if (spacer2 > 2) continue;
		    if ((abondtype & 0xf) < 1)
		     if ((dhit[nd].bondtype & 0xf) < 1)
		      e -= 0.5;
		    if (var > 7) e -= 1.0;
		    if (dl > 12) e -= 1.0;
		    if (cloop != 7) e-= 2.0;
                    if (cstem < 6) e -= 3.6;
                    else e -= 0.5; }
                 else
                  { if (cstem > 5) continue;
                    s = cpos;
                    se = cend-1;
                    while (!bp[*s][*se])
                     { s++;
                       se--; }
                    if (!stemterm[s[-1]][se[1]]) e -= 0.5;
                    e -= 0.8; }}
              else
               { if (spacer1 > 2)
		  { if (spacer1 > 3) continue;
		    if (dstem > 4) continue;
		    if (dstem < 3) continue;
		    if (tl > 15) continue;
		    if (astem < 7) e -= 1.0;
		    if (ti > 4) e -= 1.0;
		    if (cloop != 7) e-= 2.0;
		    if (tbondtype > 0x2000)
		     if (!RI[tpos[tstem-1]]) e -= 2.0;
	            e -= 1.0;
		    if (spacer2 != 1) e -= 0.5;
		    else
		     if (dhit[nd].bondtype < 0x100)
		      if (var >= 3)
                       if (var <= 5)
                        if (tstem >= 3)
			 { e += 1.0;
                           if (agcat >= 5)
                            if (wcbp[*aend1][*apos2])
                             if (!bp[aend1[-1]][*apos2])
		              if (bp[b48][dpos[dstem+1]])
                               e += 0.5; }
		    }
                 if (spacer2 > 1)
                  { if (spacer2 > 2) continue;
                    if (astem < 7)
                     if (spacer1 == 2)
                      e -= 1.0;
		    if (cloop != 7) e -= 2.0;
		    if (ea < -5.8) e -= 2.0; 
                    e -= 2.5; 
		    if (bp[b48][dpos[dstem+1]]) 
                     { if (dhit[nd].bondtype < 0x1000)
                        if (wcbp[dpos[1]][dpos[darm-2]])
                         if (wcbp[dpos[2]][dpos[darm-3]])
                          if (var < 6)
                           if (dl > 3)
                            e += 2.0; }
                    else e -= 1.0;
		     }
                 else
		  if (spacer2 < 1)
		   { if (spacer2 < 0) continue;
                     if (var > 6) continue;
		     if (dstem > 4) continue;
		     if (dhit[nd].stem_energy < -4.3) continue;
                     if (astem < 7)
                      if (spacer1 == 2)
                       e -= 1.0;
		     if (cloop != 7) e-= 2.0;
		     e -= mtBONDSTAB; }
                 if (cstem > 5)
                  if ((!gt[*cpos][cend[-1]]) || astem8) e-= mtBONDSTAB;  }

  /* very large or very small dloop */

            if (dl < 3) e -= 2.0;
	    if (dl > 14) e -= 2.0;
	    else
	     if (dl > 13)
	      { if (dhit[nd].bondtype >= 0x100) e -= 2.0;
	        else e -= 1.0; }
	     else
	      if (dl > 12)
	       { if (dhit[nd].bondtype >= 0x1000) e -= 2.0;
	         else e -= 1.0; }
              else
               if (dl > 11)
	        if (dhit[nd].bondtype >= 0x10000) e -= 2.0;

  /* tertiary interactions in darm */

             b8 = dpos[-2];
             b9 = dpos[-1];
             if (dl > 2)
               { if (dl > 5)
                  if (!stackbp[dpos[dstem+1]][b48]) e -= 1.0;
                 if (!stackbp[b8][dpos[dstem]]) e-= 0.25;
                 if (!stackbp[b8][dpos[dstem+dl-1]]) e -= 0.25; }
              if (!bp[b9][dpos[2]])
               if (!bp[b9][dpos[darm-3]])
                e -= 1.0;

  /* TR motif at b8-9 */

              if (RI[b9])
               { if (b8 == Thymine)
                  if (spacer1 == 2)
                   if (ti < 6)
                    if (((bondtype & 0xf) > 2) || (bondtype < 0x1000) ||
	                ((tbondtype < 0x100) && (tstem > 3)))
                     if ((cbondtype & 0xf) < 5)
                    if (stembp[dpos[1]][dpos[darm-2]])
                    if (var < 6)
                     if (var > 2) e += 1.5;
		     else
                      if (tstem > 3)
			   if (cloopend[-2] == Adenine)
		        e += 1.5; }
              else
               { e -= 1.0;
                 if (b9 == Thymine)
                  if (spacer1 == 2) e -= 2.0; }
              if (e > ed)
               { ed = e;
                 ndi = nd; }}

          if (ndi < 0) continue;
          energy = 100.0 + ec + ed + et;
          dl = dhit[ndi].loop;
          dstem = dhit[ndi].stem;
          darm = 2*dstem + dl;
          dpos = dhit[ndi].pos;
          dbondtype = dhit[ndi].bondtype;
          spacer1 = (int)(dpos - aend1);
          spacer2 = (int)(cpos - dpos) - darm;
          b8 = dpos[-2];

  /* tertiary structure interaction between tloop and dloop */

          if (tl >= 3)
           if (dl >= 4)
            { di = (dl < 7)?(darm-dstem-2):(darm-dstem-3);
              ti = (tl < 9)?(tstem+2):((tl < 13)?(tstem+3):(tstem+5));
              if (ggbp[dpos[di]][tpos[ti]])
               if (ggbp[dpos[di-1]][tpos[ti-1]])
                { energy += 2.0;
                  if (spacer1 != 2)
                   if (spacer2 != 2)
                    if (dstem < 4)
                     if (tl > 7)
                      if (bp[dpos[di+1]][tpos[ti+1]])
                       energy += 4.0;
                  if (ea > -2.5)
                   if (wcbp[dpos[1]][dpos[darm-2]])
                    if (wcbp[dpos[2]][dpos[darm-3]])
                     energy += 3.0; }
              if (tl > 10)
               if (dl > 10)
                energy -= 1.0; }
           else
            if (dl == 3)
             if (wcbp[dpos[dstem+1]][b48]) energy += 1.0;

  /* small darm and tarm */

         if (tloop <= 18)
          if (tarm <= 13)
           if (dl <= 8)
            if (spacer1 == 2)
             if (spacer2 == 1)
              if (abondtype < 0x1000)
               if (tbondtype < 0x100)
                if (dbondtype < 0x200)
                { et = (mtBONDSTAB - 0.5)*(double)(5 - tstem) + 
			0.1*(double)(7-tl);
                  ed = mtBONDSTAB*(double)(4 - dstem);
                  energy += (0.8*(et + ed)); }

  /* GTTC motif on tloop */

          s = tpos + tstem;
          if (tl < 5)
           if (tl < 2) energy += G[s[-1]];
           else
           { et = (G[s[-1]] + T[*s] + T[s[1]]);
             if (tl > 3)
             if (bp[*s][s[tl-1]])
              { e = (G[*s] + T[s[1]] + T[s[2]]);
                if (e > et) et = e; }
             if (tstem < 5)
              { e = (G[s[-2]] + T[s[-1]] + T[*s] + C[s[1]]);
                if (e > et) et = e; }
             energy += et; }
         else energy += (G[s[-1]] + T[*s] + T[s[1]] + C[s[2]]);

  /* long astem */

         if (astem8)
          if (bp[apos1[0]][apos2[6]])
           if (bp[apos1[1]][apos2[5]])
            if (bp[apos1[2]][apos2[4]])
             if (bp[apos1[3]][apos2[3]])
              energy += hbem[apos1[-1]][apos2[7]];

  /* false positive supression */

         if (!RI[cend[0]]) energy -= 1.0;
         if (!RI[cpos[-1]]) energy -= 1.0;
	 if (tarm < (var + 3)) energy -= 2.0;
         if (gcv < 1.5)
          if (dbondtype > 0x10000) energy -= 2.0;
	 if (tarm > 27)
          { energy -= 1.0;
            if (spacer2 != 1) energy -= 1.0; }
	 if (dstem < 3)
	  { if (var > 5) energy -= 1.0;
            if (tloop > (dloop + 8)) energy -= 0.5; }
	 if (b8 != Thymine)
          if (dl > 3)
           if (dbondtype > 0x100)
	    if ((b8 == Cytosine) ||
               (dbondtype > 0x10000))
             if (*clooppos != Cytosine)
              if (!wcbp[dpos[dstem+1]][b48])
               energy -= 1.0;


  /* high GC false positive suppression */

	 if (gcv >= 5.1)
	  { if ((abondtype & 0xf) >= 4)
             { s1 = apos1;
               s2 = apos2 + astem;
               n = 0;
               while (--s2 >= apos2)
                if (gc[*s1++][*s2])
                 { if (++n >= 4)
		    { energy -= 2.0;
                      break; }}
                else n = 0; }
	    if ((dbondtype & 0xf) >= 4) energy -= 3.0;
	    if ((cbondtype & 0xf) >= 5) energy -= 3.5;
	    if ((tbondtype & 0xf) >= tstem) energy -= 4.0; }

  /* global stem damage level */

	 tc = tstem + dstem;
	 dtbondtype = dbondtype + tbondtype;
         mabondtype = dtbondtype + cbondtype;
         bondtype = acbondtype + dtbondtype;
         if (bondtype < 0x100) energy += 0.5;
	 if ((dtbondtype & 0xf) < 1)
          { energy -= 1.0;
            if (tc >= 10) energy -= 2.0;
	    if ((bondtype & 0xf) < 3)
             if (nbase > 75) energy -= 1.0; }
         i = (int)((bondtype >> 16) & 0xf);
         j = (int)((bondtype >> 12) & 0xf) + i;
         k = (int)((bondtype >> 8) & 0xf) + j;
	 ti = (tc > 6)?5:((tc > 5)?4:3);
         if (k > ti)
          { ev = (double)(k - ti);
            energy -= 0.5*ev;
            if (cbondtype > 0x10000)
             if (tstem < 5)
              energy -= ev;
            if (i > 0)
             if (k > 8)
              energy -= 1.5*(double)(k - 8); }

  /* low GC false positive supression */

         if (gcv < 3.5)
          if ((bondtype & 0xf) < 2)
           { if ((bondtype & 0xf) < 1) energy -= 1.0;
              if (dl > 3)
               if (var > 2)
                if (!wcbp[dpos[dstem+1]][b48]) energy -= 1.0; }

  /* small variable loop */

         if (var < 3)
	  { if (dloop > 18)
	     { if (dloop > (tloop + 2)) energy -= 1.0;
	       if (tloop > 20)
                if ((((dtbondtype >> 4) + dtbondtype) & 0xf) < 6)
	         energy -= 2.0; }
	    if (astem < 7) 
             { energy -= 1.0;
               if (agcat >= 5)
	        if (bondtype < 0x300)
	         if (gcv > 1.2)
                  if (gcv < 5.0)
                   energy += 2.0; }}
	 else

  /* NNNNNAA cloop */

          if (cloopend[-2] == Adenine)
           if (cloopend[-1] == Adenine)
            if (spacer1 > 1)
            if ((dbondtype < 0x2000) || (dloop > mt_DRLmaxlength))
             { if (abondtype < 0x100) energy += 1.0;
               else
		if (cbondtype < 0x100) energy += 1.0;
	        else
                 if (tstem >= 5)
		  if (tbondtype < 0x100)
                   { energy += 1.0;
                     if (*clooppos == Cytosine)
		      if (clooppos[1] == Thymine)
                       if (dbondtype < 0x100)
		       energy += 0.5;

                     if (cgcat >= 3)
                     if ((tbondtype & 0xf) > 0)
                     if (ggbp[dpos[dstem+1]][b48])
                     if (wcbp[dpos[1]][dpos[darm-2]])
                     if (tl < 10)
                     if (spacer1 == 2)
	              if (spacer2 == 1)
                       if (dl > 2)
                        if (var >= 2)
                         if (var < 6)
                          { if (agcat >= 6) energy += 3.0;
			    else
			     if (agcat >= 5)
                              if (cgcat >= 4)
                               if (dbondtype < 0x100)
                                if (*s == Thymine)
                                 if (s[-1] == Guanine)
                                  if (s[1] == Thymine)
                                   energy += 3.0; }}}

  /* large tloop */

         if (tl > 12)
	  { if (tbondtype > 0x10000) energy -= 2.0;
            if (agcat < 5)
             if (spacer1 != 2)
              if (spacer2 != 1)
               energy -= 1.0; }

  /* find exceptions */

	 if (energy < dtthresh)
          { if (!mtxdetect) continue;
            if (energy < (thresh - 12.0)) continue;
            if (energy < (dtthresh - 12.0)) continue;
	    if (nbase > 75) continue;
	    if (dstem > 4) continue;
	    if (dstem < 3) continue;
	    if (astem < 7)
	     if (acbondtype > 0x21000)
              continue;
            if (var > 5)
             { if (var > 6) continue;
               if (tarm < 12) continue; }
	    if (gcv <= 1.2)
             { if (gcv < 0.9) continue;
               if ((mabondtype & 0xf) < 1) continue; }
	    if (tl > 9)
             { if (tl > 13) continue;
               if (!wcbp[dpos[1]][dpos[darm-2]])
                continue; }
	    if (dl > 7)
             { if (bondtype > 0x20000)
                if (dloop > (tloop + 4)) continue;
	       if (dl > 10)
                { if (dl > 12)
                   if (abondtype > 0x1000)
		    continue;
	          if (tbondtype > 0x200) continue;
                  if (tt[*tpos][apos2[-1]]) continue;
                  if (var > 5) continue;
	          if (dloop > (tloop + 8))
		    if (bondtype > 0x10000) continue;
		  if (astem < 7) continue; }}
	    if (RI[clooppos[1]]) continue; 
	    b9 = dpos[-1];
	    if (cstem >= 6)
	     { if (cbondtype > 0x200) continue;
               if (var < 3) continue; 
	       if (YI[b9]) continue; }
	    if (cloop != 7) continue;
            if (ds == MAMMAL_MT) continue;
            if (mabondtype < 0x400)
             { if ((b8 == Thymine) || (mabondtype < 0x300))
		if (ea < -5.45)
	         if (chit[nc].stem_energy > -3.2)
		  if (dbondtype < 0x200)
                  if (spacer1 > 1)
                   if ((spacer2 == 1) || (mabondtype < 0x100))
                    if ((spacer1 < 3) || (tstem > 3) ||
			(tbondtype < 0x100))
                     if ((spacer1 < 3) ||
		      ((var > 2) && (var < 6) && (tbondtype < 0x2000) &&
		       (tl < 10)))
                    if (dstem < 5)
                   if (var >= 2)
 	            if (dl > 2)
                     if (tl < 15)
                      if ((b8 != Cytosine) ||
		        (*clooppos == Cytosine))
                     if (RI[b9])
                      if (*clooppos != Adenine)
                       if (clooppos[1] == Thymine)
			if (RI[cloopend[-2]])
                         { s1 = apos1;
                           s2 = apos2 + astem;
                           n = 0;
                           while (--s2 >= apos2)
                            if (wcbp[*s1++][*s2])
                             { if (++n >= 3) break; }
                            else n = 0;
	                   if (n >= 3)
			    { energy += 3.0;
                              if ((abondtype & 0xf) > 0) energy += 2.0;
                              if (bp[dpos[dstem+1]][b48])
                                if (wcbp[dpos[1]][dpos[darm-2]])
                                 if (var <= 5)
                                  energy += 1.0; }
			   if (dtbondtype < 0x200)
                           if (agcat < 2)
                            if (wcbp[dpos[dstem+1]][b48])
                             if (wcbp[dpos[1]][dpos[darm-2]])
                              if (wcbp[dpos[2]][dpos[darm-3]])
			      if (gcv > 1.2)
                              if (var <= 5)
			       if (tstem >= 3)
			        if (dstem >= 3)
                                 if (tl > 3)
                                  if (tl < 9)
                                   if (dl < 9)
                                    if (spacer1 == 2)
				     energy += 10.0; }
               if ((tbondtype & 0xf) > 0)
               if (mabondtype < 0x300)
	        { if (mabondtype < 0x100)
                   { if ((spacer1 < 3) || (tstem > 2))
                    if (var > 0)
                     if (YI[*clooppos])
                      if ((spacer2 > 0) || (clooppos[1] == Thymine))
	               energy += 2.5; }
	        else
                 if ((dbondtype & 0xf) > 0)
	          if (b9 != Cytosine)
                   if (var <= 7)
		    if (spacer2 == 1)
                     if (tarm < 22)
                      if (gcv > 1.2)
                     if (dstem >= 4)
		      { if (tstem >= 5)
                       energy += 5.0;
	              else
	               if (tstem >= 3)
	                if (tbondtype < 0x100)
	                 energy += 1.0; }
	             else
		      if (tstem >= 5)
                        energy += 1.0; }
             else
              if ((dbondtype & 0xf) > 0)
              { if (tstem >= 5)
		 if (s[-1] == Guanine)
		  if (*s == Thymine)
		   if (s[1] == Thymine)
                    if (*clooppos == Cytosine)
                     if (clooppos[1] == Thymine)
			  if (cloopend[-2] == Adenine)
                       if (cloopend[-1] == Adenine)

			energy += 1.0;
	         if (bondtype < 0x1000)
	          if (cbondtype < 0x100)
	           energy += 1.0; }}
			
	    if (tstem >= 5)
		 if (*clooppos == Cytosine)
		 {
             if (dl > 3)
             if (dtbondtype < 0x200)
              if ((tbondtype & 0xf) > 0)
              if (clooppos[1] == Thymine)
               {
	        if (clooppos[2] == Thymine)
		  if (clooppos[3] == Adenine)
              if (clooppos[4] == Cytosine)
              if (clooppos[5] == Adenine)
			  if (cloop == 7)
		      energy += 0.5;

                if (cgcat >= 4)
		 if (wcbp[dpos[1]][dpos[darm-2]])
                  if (bp[dpos[dstem+1]][b48])
                     if (tl < 10)
                      if (var < 6)
	               if (spacer1 == 2)
	                if (spacer2 == 1)
		        if (dstem >= 3)
		         energy += 3.0;

	       }
			
			
            if (clooppos[1] == Cytosine)
              if (clooppos[2] == Cytosine)
		  if (clooppos[3] == Adenine)
	        if (clooppos[4] == Thymine)
  		  if (s[-1] == Guanine)
		   if (*s == Thymine)
		    if (s[1] == Thymine)
		     energy += 1.0;
		}	

		
		if (RI[b9])
		{	
		if (b8 == Thymine)
		{
		if (clooppos[1] == Thymine)
                {
		if (cloopend[-2] == Adenine)
		{
		
		if (wcbp[dpos[1]][dpos[darm-2]])
		{
                if (*clooppos == Cytosine)
	        {
	         if (abondtype < 0x200)
		   energy += 1.0;

                 if (bondtype < 0x10000)
                  if (dtbondtype < 0x200)
                   if (agcat >= 3)
                    if (cgcat >= 4)
                     if (tl < 10)
                      if (var < 6)
	               if (spacer1 == 2)
	                if (spacer2 == 1)
                          if (tstem >= 3)
		            energy += 3.0;
		}

			
                if (tstem >= 5)
		 if (s[-1] == Guanine)
		  if (*s == Thymine)
		   if (s[1] == Thymine)
                    { energy += 1.0;
                      if (tl >= 5)
	               if (spacer1 == 2)
	                if (spacer2 == 1)
                         if (tbondtype < 0x100)
                          if (wcbp[dpos[dstem+1]][b48])
                           energy += 3.0; }

                if (tstem >= 3)
                if (tl < 10)
	        if (spacer1 == 2)
	        if (spacer2 == 1)
		 if (RI[cloopend[-1]])
                  if (dl > 2)
                   if (var >= 2)
                  if (var < 6)
                  if (ggbp[dpos[dstem+1]][b48])
		  {
                   if (dtbondtype < 0x100)
                    { energy += 3.5;
                      if ((bondtype & 0xf00) == 0)
                       if (*clooppos == Cytosine)
                         energy += 1.5; }


                   if (bondtype < 0x10000)
	            if (tstem > 2)
                     if (tbondtype < 0x200) energy += 2.5;

		   if (abondtype < 0x100)
		    if (wcbp[dpos[2]][dpos[darm-3]])
                     energy += 3.0;

                   if (tbondtype < 0x100)
                   if (agcat >= 6)
                   if (tstem >= 5)
                   if ((tbondtype & 0xf) > 0)
                   if (RI[cloopend[-1]])
                   if (cgcat >= 4)
                    energy += 2.0;
		  }
		  else
                   if (!ggstembp[*tpos][apos2[-1]])
                    if (wcbp[dpos[dstem+1]][*tpos])
                     energy += 1.5;
        }
			
	    if ((abondtype & 0xf) < 1)
	    if (abondtype < 0x100)
             if (gcv > 1.2)
	      if (dl > 3)
                if (bp[dpos[dstem+1]][b48])
                if (spacer1 == 2)
                 if (spacer2 == 1)
             if (*clooppos == Cytosine)
		      energy += 5.0;

	    if (cbondtype < 0x100)
             if (tbondtype < 0x100)
              if (tstem >= 3)
	      if (dl > 3)
               if (var < 6)
                if (bp[dpos[dstem+1]][b48])
                if (spacer1 == 2)
                 if (spacer2 == 1)
		  energy += 2.5;

        }

        if (stembp[dpos[dstem+1]][b48])
		{

             if (*clooppos == Thymine)
                if (cloopend[-2] == Guanine)
              if (clooppos[2] == Guanine)
	       if (clooppos[3] == Thymine)
	        if (clooppos[4] == Guanine)
		  if (dl > 2)
              energy += 1.0;
				
	    if (cbondtype < 0x100)
             if (dbondtype < 0x10000)
                  if (wcbp[dpos[1]][dpos[darm-2]])
		    if (var < 6)
		     if (tstem >= 3)
			if (gcv >= 1.2)
			  if (dl > 3)
                energy += 1.0;

	    if (tstem >= 5)
             if (dtbondtype < 0x200)
              if (*clooppos == Cytosine)
                if (spacer1 == 2)
                 if (spacer2 == 1)
		 if (RI[cloopend[-2]])
                  energy += 0.5;
		}	

		}		
				

	     if (tstem > 2)
              if (tarm < 28)
                if (spacer1 == 2)
                 if (spacer2 == 1)
                  if (dl > 3)
                   if (j < 1)
		    if (k > ti)
                   if (ggstembp[dpos[dstem+1]][b48])
			energy += 2.5;
			
            if (dtbondtype < 0x100)
              if ((tbondtype & 0xf) > 0)
               if (bp[dpos[dstem+1]][b48])
                 if (b9 == Adenine)
                  if ((dbondtype & 0xf) > 0) energy += 2.0;
	          else
                   if (spacer2 == 1)
	            energy += 0.5;
	    if (cloopend[-2] == Adenine)
	     if (cloopend[-1] == Adenine)
              if (cbondtype < 0x2000)
               if (spacer1 > 1)
 	       if (dl > 2)
                if (var < 6)
                  energy += 0.75;
				
		
		}
		
            if (var > 2)
             if (dl > 2)
              {
              if (cbondtype < 0x200)
               if (((mabondtype & 0xf) > 3) || (bondtype < 0x1000))
		{
                if (bp[dpos[dstem+1]][b48])
                 energy += 1.0;

                if (cbondtype < 0x100)
		if (dbondtype < 0x100)
	        if (bp[b8][dpos[dstem]])
                if (bp[b8][dpos[darm-dstem-1]])
                if (var < 6)
                if (tstem >= 3)
                if (tl < 10)
	        if (spacer1 == 2)
	        if (spacer2 == 1)
                if (clooppos[1] == Thymine)
		 if (cloopend[-2] == Adenine)
                  energy += 3.0;
		}

              if (clooppos[1] == Thymine)
		   if (RI[cloopend[-2]])
		{
                if (*clooppos == Cytosine)
		{

	      if (dtbondtype < 0x200)
                if (agcat >= 3)
                if (cgcat >= 4)
                if (var < 6)
                if (tstem >= 3)
                if (tl < 10)
	        if (spacer1 == 2)
	        if (spacer2 == 1)
                 { if (abondtype > 0x20000)
                    if (bp[dpos[dstem+1]][b48])
		     energy += 7.0;
                   if (agcat >= 6) energy += 2.0; }
	
	        if ((bondtype & 0xf00) == 0)
		 if (gcv > 5.0)
  		  if (s[-1] == Guanine)
		   if (*s == Thymine)
                    if (tstem >= 5)
                    if (var < 6)
                    if (tl < 10)
	            if (spacer1 == 2)
	            if (spacer2 == 1)
		     energy += 2.0;

		if (abondtype < 0x100)
                 if (cbondtype < 0x10000)
                  if (bp[dpos[dstem+1]][b48])
                   if (cgcat >= 4)
                    if (tstem >= 3)
                     if (var < 6)
                      if (tstem >= 3)
                       if (tl < 10)
	                if (spacer1 == 2)
	                 if (spacer2 == 1)
			  energy += 1.5;
	        }


	      if (dtbondtype < 0x100)
               if (agcat >= 4)
               if (cgcat >= 4)
                if (var < 6)
                if (tstem >= 3)
                if (tl < 10)
                 { if (spacer1 == 2)
	            { if (abondtype < 0x3000)
                       if (stackbp[dpos[dstem+1]][b48])
		        energy += 3.0;
		      if (b8 == Thymine)
  		       if (s[-1] == Guanine)
		        if (*s == Thymine)
		         if (s[1] == Thymine)
		          energy += 3.5; }
                   if (agcat >= 6)
		    if (YI[*clooppos])
  		     if (s[-1] == Guanine)
		      if (*s == Thymine)
                       if ((dtbondtype & 0xf) > 0)
                        energy += 3.0; }

              if (mabondtype < 0x10000)
               if (dtbondtype < 0x400)
                if (agcat >= 5)
                 if (cgcat >= 3)
                  if (tl < 10)
                   if (var < 6)
	            if (spacer1 == 2)
	             if (spacer2 == 1)
		     {
		
                     if (dtbondtype < 0x200)
                      if (cbondtype < 0x300)
                      if (bondtype < 0x10000)
                       if (tstem >= 3)
		         energy += 1.0;

		     if (tstem >= 5)
                      if (s[-1] == Guanine)
		       energy += 4.0;
		     }
		}
	      }
		
		}
		else
                 if (bondtype < 0x10000)
                 if (mabondtype < 0x500)
                 if (dbondtype < 0x100)
		 if (b8 == Thymine)
                 if (agcat >= 4)
		 if (clooppos[1] == Thymine)
		 if (cloopend[-2] == Adenine)
		 if (cloopend[-1] == Adenine)
	         if (spacer1 == 2)
	         if (spacer2 == 1)
                 if (dstem >= 3)
                 if (tstem >= 5)
                 if (dl > 2)
                 if (tl < 10)
                 if (var < 6)
		  energy += 7.0;

           if (agcat >= 5)
            if (cgcat >= 4)
             if ((acbondtype & 0xf) >= 3)
	     { if (tbondtype < 0x100)
                if ((dbondtype & 0xf) > 0)
                if ((((dbondtype >> 4) + dbondtype) & 0xf) >= 3)
                if (wcbp[dpos[dstem+1]][b48])
                if (b8 == Thymine)
                  if (RI[b9])
		   if (clooppos[1] == Thymine)
                   if (YI[*clooppos])
		   if (RI[cloopend[-2]])
		   if (RI[cloopend[-1]])
	           if (spacer1 == 2)
	           if (spacer2 == 1)
                   if (dl > 2)
                   if (tl < 10)
                   if (var < 6)
                   if (var > 2)
                    energy += 6.0;
               if (cgcat >= 5)
                if (abondtype < 0x10000)
                 if (bp[dpos[dstem+1]][b48])
		   if (clooppos[1] == Thymine)
                   if (YI[*clooppos])
		   if (RI[cloopend[-2]])
                   if (dl > 2)
                   if (tl < 10)
                   if (var < 6)
                   if (var > 2)
                    energy += 6.0; }

	    if (energy >= dtthresh)
             energy -= (0.9*(energy - dtthresh) + 5.0);
	    else continue; }

  /* remember fully formed mttRNA gene if threshold reached */

            if (energy < thresh) continue; 
	    te.energy = energy;
            thresh = energy;
            te.ps = apos1;
            te.spacer1 = spacer1;
         te.dstem = dstem;
         te.dloop = dl;
            te.spacer2 = spacer2;
            te.cstem = cstem;
         te.cloop = cloop;
	    te.var = var;
	        te.varbp = (var > 17)?varbp:0;
            te.tstem = tstem;
            te.tloop = tl;
         k = astem + spacer1 + darm + spacer2;
         te.anticodon = k + cstem + 2;
         te.nintron = 0;
         te.intron = 0;
         te.nbase = k + carm + var + 2*tstem + tl;
	 tastem = astem;
	 tastem8 = astem8;
	 tastem8d = astem8d;
         } }

  /* for highest energy mttRNA gene */ 
  /* decide astem length, look for NCCA acceptor tail */
  /* and calculate total length */

    if (te.ps)
     {    apos2 = te.ps + te.nbase;
          if (extastem)
           if (tastem8d)
            { te.astem1 = 8;
              te.astem2 = 8;
              te.ps--;
	      te.nbase++;
	      te.anticodon++;
              as = aatail(apos2+8,&aext,sw); }
           else
            { te.astem1 = tastem;
              te.astem2 = tastem;
              as = aatail(apos2+tastem,&aext,sw);
              if (tastem8)
               { as8 = aatail(apos2+8,&aext8,sw);
                 if (as8 >= as)
                  { te.ps--;
	            te.nbase++;
	            te.anticodon++;
                    te.astem1 = 8;
                    te.astem2 = 8;
                    as = as8;
                    aext = aext8; }}}
          else
           { te.astem1 = tastem;
             te.astem2 = tastem;
             as = aatail(apos2+tastem,&aext,sw); }
          if (as < 2) aext = 1;
          te.nbase += te.astem2;
          nbasefext = te.nbase + ASTEM2_EXT;
          te.nbase += aext;

  /* store mttRNA gene if there are no */
  /* higher energy overlapping mttRNA genes */

          te.start = (long)(te.ps - seq);
          if (tn = find_slot(d,&te,&nto,sw))
           { base_copy3(te.ps,te.seq,nbasefext);
             base_copy3(te.ps,te.eseq,nbasefext);
             te.aatail = aext;
             *tn = te;  }}
     }
  return(nto); }



void init_gene(gene *t, int nstart, int nstop)
{ int i;
  for (i = nstart; i < nstop; i++)
   { t[i].energy = -1.0;
     t[i].genetype = noGENE;
     t[i].tps = 0;
     *(t[i].name) = '\0'; }}

void disp_ftable_entry(FILE *f, int n[], int i, int m, csw *sw)
 { if (m > 0)
            switch(sw->geneticcode)
              { case METAZOAN_MT:
                        if (i < 2) fprintf(f," %-4s %-4d",aa(n,sw),m);
                        else fprintf(f," %-18s %-4d",aa(n,sw),m);
                        break;
                case STANDARD:
                case VERTEBRATE_MT:
                default:
                        fprintf(f," %-4s %-5d",aa(n,sw),m);
                        break; }
           else
            switch(sw->geneticcode)
              { case METAZOAN_MT:
                        if (i < 2) fprintf(f," %-4s     ",aa(n,sw));
                        else fprintf(f," %-18s     ",aa(n,sw));
                        break;
                case STANDARD:
                case VERTEBRATE_MT:
                default:
                        fprintf(f," %-4s      ",aa(n,sw));
                        break; }}


void disp_freq_table(gene *t, int nt, csw *sw)
{ int i,j,k,m,ambig,*s,c1,c2,c3,c[3],n[3],table[4][4][4];
  static int cgflip[4] = { 0,2,1,3 };
  FILE *f = sw->f;
  ambig = 0;
  for (i = 0; i < 4; i++)
   for (j = 0; j < 4; j++)
    for (k = 0; k < 4; k++)
     table[i][j][k] = 0;
  for (i = 0; i < nt; i++)
   if (t[i].energy >= 0.0)
    if (t[i].genetype == tRNA)
     if (t[i].cloop == 7)
     { s = t[i].seq + t[i].anticodon;
       c1 = *s;
       c2 = s[1];
       c3 = s[2];
       if ((c1 >= Adenine) && (c1 <= Thymine))
        if ((c2 >= Adenine) && (c2 <= Thymine))
         if ((c3 >= Adenine) && (c3 <= Thymine))
          table[*s][s[1]][s[2]]++;
         else ambig++;
        else ambig++;
       else ambig++; }
     else ambig++;
  fprintf(f,"tRNA Anticodon Frequency\n");
  for (j = 0; j < 4; j++)
   { n[2] = cgflip[j];
     for (k = 0; k < 4; k++)
      { n[1] = cgflip[k];
        for (i = 0; i < 4; i++)
         { n[0] = cgflip[i];
           fprintf(f,"%c%c%c",cpbase(n[0]),cpbase(n[1]),cpbase(n[2]));
           m = table[n[0]][n[1]][n[2]];
          disp_ftable_entry(f,n,i,m,sw); }
        fputc('\n',f); }}
  if (ambig > 0) fprintf(f,"Ambiguous: %d\n",ambig);
  fprintf(f,"\ntRNA Codon Frequency\n");
  for (i = 0; i < 4; i++)
   { n[0] = 3 - cgflip[i];
     for (j = 0; j < 4; j++)
      { n[1] = 3 - cgflip[j];
        for (k = 0; k < 4; k++)
         { n[2] = 3 - cgflip[k];
           fprintf(f,"%c%c%c",cpbase(n[0]),cpbase(n[1]),cpbase(n[2]));
           c[0] = 3 - n[2];
           c[1] = 3 - n[1];
           c[2] = 3 - n[0];
           m = table[c[0]][c[1]][c[2]];
           disp_ftable_entry(f,c,k,m,sw); }
        fputc('\n',f); }}
  if (ambig > 0) fprintf(f,"Ambiguous: %d\n",ambig);
  fputc('\n',f); }

void disp_energy_stats(data_set *d, gene t[], int nt, csw *sw)
{ int i,n[NS],genetype,introns,nintron,trna,mtrna,ntv,nd,nps;
  double gc,gcmin[NS],gcmax[NS];
  static char genetype_name[NS][30] =
   { "tRNA genes","","","","Overall" };
  FILE *f = sw->f;
  mtrna = sw->mtrna;
  trna = sw->trna | mtrna;
  nps = 0;
  if (mtrna)
   { ntv = 0;
     nd = 0; }
   introns = 0;
  for (i = 0; i < NS; i++)
   { n[i] = 0;
     gcmin[i] = 1.0;
     gcmax[i] = 0.0; }
  for (i = 0; i < nt; i++)
   if (t[i].energy >= 0.0)
    { n[NS-1]++;
      genetype = t[i].genetype;
      n[genetype]++;
      if (ts[i].energy < 100.0) nps++;
      if (genetype == tRNA)
       { if (mtrna)
   { if (t[i].tstem == 0) ntv++;
     if (t[i].dstem == 0) nd++; }
         gc = gc_content(t+i);
         if (gc < gcmin[genetype]) gcmin[genetype] = gc;
         if (gc > gcmax[genetype]) gcmax[genetype] = gc; }}
  fputc('\n',f);
  fputc('\n',f);
  if (sw->repeatsn)
   if ((n[0] + n[1]) > 0)
     fprintf(f,"%s\n\n",d->seqname);
  sw->ngene[0] += n[0];
  if (n[0] > 3) disp_freq_table(t,nt,sw);
  if (n[0] > 1) 
   { fprintf(f,"Number of %s = %d\n",genetype_name[0],n[0]);
     if (sw->tvloop)
      fprintf(f,"Number of TV replacement loop tRNA genes = %d\n",
             ntv);
     fprintf(f,"Number of D replacement loop tRNA genes = %d\n",
     nd); 
     if (n[0] > 1)
      fprintf(f,"tRNA GC range = %2.1f%% to %2.1f%%\n",
                 gcmin[0]*100.0,gcmax[0]*100.0); }
  for (i = 2; i < NS-1; i++)
   if (n[i] > 1)
    { sw->ngene[i] += n[i];
      fprintf(f,"Number of %s = %d\n",genetype_name[i],n[i]); }
  if (sw->reportpseudogenes)
   if (nps > 0) 
    if (n[NS-1] > 1)
     fprintf(f,"Number of possible pseudogenes = %d\n",nps);
  fputc('\n',f);
  fputc('\n',f); }

int gene_sort(data_set *d, gene *ts, int nt, int sort[])
{ int i,n,j,k;
  long starti,startj,stopi,stopj;
  n = 0;
  for (i = 0; i < nt; i++)
   if (ts[i].energy >= 0.0)
    sort[n++] = i;
  i = -1;
  while (++i < (n-1))
   { j = i;
     while (++j < n)
      { starti = ts[sort[i]].start;
 startj = ts[sort[j]].start;
      if (starti > startj)
       { k = sort[i];
  sort[i] = sort[j];
  sort[j] = k; }
      else
       if (starti == startj)
 { stopi = ts[sort[i]].stop;
          stopj = ts[sort[j]].stop;
   if (stopi < starti) stopi += d->psmax;
   if (stopj < startj) stopj += d->psmax;
          if (stopi < stopj)
           { k = sort[i];
      sort[i] = sort[j];
      sort[j] = k; }}}}
   return(n); }



int iamatch(gene *t, csw *sw)
{ char key[5],*k,s[100];
  if (!(k = softstrpos(t->name,"TRNA-"))) return(-1);
  copy3cr(k+5,key,3);
  name(t,s,1,sw);
  if (softstrpos(s,key)) return(1);
  return(0); }

void disp_gene_set(data_set *d, gene *ts, int nt, csw *sw)
{ int i,j,n,vsort[NT],*sort;
  char m[MATX][MATY],s[20];
  gene *t;
  FILE *f = sw->f;
  if (nt <= NT)
   sort = vsort;
  else
   { sort = (int *)malloc(nt*sizeof(int));
     if (sort == NULL)
      { fprintf(stderr,"Not enough memory to sort genes\n");
        exit(1); }}
  n = gene_sort(d,ts,nt,sort);
  j = sw->trna_struct[49];
  for (i = 50; i <= 55; i++) j += sw->trna_struct[i];
  if (j != ((sw->trna_struct[0] << 4) + 9)) return;
  if (sw->libflag < 2)
   { if (n > 0)
      for (j = 0; j < n;)
       { i = sort[j++];
         t = ts + i;
         t->energy = nenergy(t,sw);
         switch(t->genetype)
          { case tRNA:
                    init_matrix(m);
                    disp_gene(t,m,sw);
                    sprintf(s,"%d.",j);
                    xcopy(m,0,32,s,length(s));
                    disp_matrix(f,m,MATY);
	            if (sw->matchacceptor)
                     if (iamatch(t,sw) == 0)
		      { fprintf(f,"Iso-acceptor mismatch\n");
			sw->iamismatch++; }
                    overlap(d,ts,sort,n,i,sw);
                    if (sw->seqdisp) disp_seq(f,t,sw);
                    break; }}
     else
      if (*(d->seqname) != '\0')
       fprintf(f,"\nNothing found in %s\n\n\n",d->seqname);
      else
       fprintf(f,"\nNothing found\n\n\n"); }
  if (sw->libflag > 0)
   { for (i = 0; i < n; i++)
      write_to_library(f,ts + sort[i],sw); }
  disp_energy_stats(d,ts,nt,sw);
  if (nt > NT) free((void *)sort); }


void batch_gene_set(data_set *d, gene *ts, int nt, csw *sw)
{ int i,j,n,vsort[NT],nspaces,caps,*sort;
  gene *t;
  FILE *f = sw->f;
  if (nt <= NT)
   sort = vsort;
  else
  {  sort = (int *)malloc(nt*sizeof(int));
     if (sort == NULL)
      { fprintf(stderr,"Not enough memory to sort genes\n");
        exit(1); }}
  n = gene_sort(d,ts,nt,sort);
  j = sw->trna_struct[49];
  for (i = 50; i <= 55; i++) j += sw->trna_struct[i];
  if (j != ((sw->trna_struct[0] << 4) + 9)) return;
  if (sw->libflag < 2)
   if (sw->batch >= 2)
    { nspaces = (sw->batch & 0x4);
      caps = (sw->batch & 0x10);
      if (sw->batch & 0x8) 
       for (i = 0; i < n; i++)
        disp_fasta_seq(f,ts + sort[i],d->ns+1,i+1,nspaces,caps,sw);
      else
       for (i = 0; i < n; i++)
        disp_fasta_seq(f,ts + sort[i],0,0,nspaces,caps,sw); }
   else  
   { fprintf(f,"%d genes found\n",n);
     for (j = 0; j < n; j++)
      { fprintf(f,"%-3d ",j+1);
        t = ts + sort[j];
        t->energy = nenergy(t,sw);
        switch(t->genetype)
         { case tRNA:  disp_batch_trna(f,t,sw);
                       break;
           default:    break; }}}
  if (sw->libflag > 0)
   { for (i = 0; i < n; i++)
      write_to_library(f,ts + sort[i],sw); }
  if (nt > NT) free((void *)sort); }




void iopt_fastafile(data_set *d, csw *sw)
{ int i,nt,flag,len,anticodon;
  int *s,*sf,*se,*sc,*swrap;
  int seq[2*LSEQ+WRAP+1],cseq[2*LSEQ+WRAP+1],wseq[2*WRAP+1];
  long gap,start,rewind,drewind,psmax,tmaxlen,vstart,vstop;
  char c1,c2,c3;
  static char trnatypename[3][25] =
   { "Metazoan mitochondrial","Cytosolic","Mammalian mitochondrial" };
  static char genecodename[NGENECODE][50] =
   { "composite Metazoan Mitochondrial",
     "standard",
     "Vertebrate Mitochondrial",
     "Yeast Mitochondrial",
     "Mold/Protozoan/Coelenterate Mitochondrial",
     "Invertebrate Mitochondrial",
     "Ciliate",
     "deleted -> standard",
     "deleted -> standard",
     "Echinoderm/Flatworm Mitochondrial",
     "Euplotid",
     "Bacterial/Plant Chloroplast",
     "Alternative Yeast",
     "Ascidian Mitochondrial",
     "Alternative Flatworm Mitochondrial",
     "Blepharisma",
     "Chlorophycean Mitochondrial",
     "deleted -> standard",
     "deleted -> standard",
     "deleted -> standard",
     "deleted -> standard",
     "Trematode Mitochondrial",
     "Scenedesmus obliquus Mitochondrial",
     "Thraustochytrium Mitochondrial" };
  FILE *f = sw->f;
  ts = (gene *)malloc(NT*sizeof(gene));
  if (ts == NULL)
   { fprintf(stderr,"Not enough memory available to store detected genes\n");
     exit(1); }
  sw->space = NT;
  init_trna(f,sw);
  fprintf(f,"\nPlease reference the following paper if you use this\n");
  fprintf(f,"program as part of any published research.\n\n");
  fprintf(f,"Laslett, D. and Canback, B. (2008) ARWEN: a\n");
  fprintf(f,"program to detect tRNA genes in metazoan mitochondrial\n");
  fprintf(f,"nucleotide sequences\n");
  fprintf(f,"Bioinformatics, 24(2); 172-175.\n\n\n");
  fprintf(f,"Searching for %s tRNA genes\n",trnatypename[sw->discrim]);
  if (!sw->tvloop)
    fprintf(f,"TV replacement loop tRNA genes not detected\n");
  if (sw->linear)
   fprintf(f,"Assuming linear topology, search will not wrap around ends\n");
  else
   fprintf(f,"Assuming circular topology, search wraps around ends\n");
  if (sw->both == 2)
   fprintf(f,"Searching both strands\n");
  else
   if (sw->both == 1)
    fprintf(f,"Searching complementary (antisense) strand only\n");
   else
    fprintf(f,"Searching single (sense) strand only\n");
  fprintf(f,"Using %s genetic code\n",genecodename[sw->geneticcode]);
  if (sw->ngcmod > 0)
      { fprintf(f,"Specified modifications:\n");
        for (i = 0; i < sw->ngcmod; i++)
         { anticodon = sw->gcmod[i];
           c1 = cpbase(Thymine - (anticodon & 0x3));
           c2 = cpbase(Thymine - ((anticodon >> 2) & 0x3));
           c3 = cpbase(Thymine - ((anticodon >> 4) & 0x3));
           fprintf(f,"%c%c%c = %s\n",c1,c2,c3,
                   aaname[aamap[sw->geneticcode][anticodon]]); }}
  fputc('\n',f);
  fputc('\n',f);
  rewind = MAXTAGDIST + 20;
  if (sw->trna | sw->mtrna)
   { tmaxlen = MAXTRNALEN + sw->maxintronlen;
     if (rewind < tmaxlen) rewind = tmaxlen; }
  sw->loffset = rewind;
  sw->roffset = rewind;
  drewind = 2*rewind;
  d->ns = 0;
  d->nextseq = 0L;
  while (d->nextseq >= 0L)
   { d->seqstart = d->nextseq;
     if (!seq_init(d)) break;
     psmax = d->psmax;
     if (sw->verbose)
      { fprintf(stderr,"%s\n",d->seqname);
        fprintf(stderr,"%ld nucleotides in sequence\n",psmax);
        fprintf(stderr,"Mean G+C content = %2.1f%%\n",100.0*d->gc);
        fprintf(stderr,"Using %s genetic code\n",genecodename[sw->geneticcode]);
        if (sw->ngcmod > 0)
            { fprintf(stderr,"Specified modifications:\n");
              for (i = 0; i < sw->ngcmod; i++)
               { anticodon = sw->gcmod[i];
                 c1 = cpbase(Thymine - (anticodon & 0x3));
                 c2 = cpbase(Thymine - ((anticodon >> 2) & 0x3));
                 c3 = cpbase(Thymine - ((anticodon >> 4) & 0x3));
                 fprintf(stderr,"%c%c%c = %s\n",c1,c2,c3,
                         aaname[aamap[sw->geneticcode][anticodon]]); }}}
     fprintf(f,"%s\n",d->seqname);
     fprintf(f,"%ld nucleotides in sequence\n",psmax);
     fprintf(f,"Mean G+C content = %2.1f%%\n",100.0*d->gc);
     init_gene(ts,0,NT);
     nt = 0;
     flag = 0;
     start = 1L;
     se = seq;
     if (sw->linear)
      { for (i = 0; i < rewind; i++) *se++ = NOBASE;
        start -= rewind; }
     else
      { if (psmax <= drewind)
         { gap = drewind - psmax;
           sc = se + gap;
           while (se < sc) *se++ = NOBASE;
           swrap = wseq;
           sc = se + psmax;
           while (se < sc)
            { *se = move_forward(d);
              *swrap++ = *se++; }
           sc = swrap + gap;
           while (swrap < sc) *swrap++ = NOBASE;
                   swrap = wseq;
                   sc = swrap + psmax;
                   while (swrap < sc) *se++ = *swrap++;
                   swrap = wseq;
                   sc = swrap + drewind;
                   while (swrap < sc) *se++ = *swrap++;
                   sw->loffset = drewind;
                   sw->roffset = drewind;
                   start -= drewind;
                   flag = 1;
                   goto SH; }
        else
         { swrap = wseq;
           sc = seq + drewind;
           while (se < sc)
            { *se = move_forward(d);
              *swrap++ = *se++; }}}
     sc = seq + LSEQ;
     NX:
     while (se < sc)
      { if (d->ps >= psmax)
         { if (sw->linear)
            for (i = 0; i < rewind; i++) *se++ = NOBASE;
           else
            { sc = wseq + drewind;
              swrap = wseq;
              while (swrap < sc) *se++ = *swrap++; }
           flag = 1;
           break; }
                else *se++ = move_forward(d); }
         SH:
     len = (int)(se - seq);
     if (sw->verbose)
      { vstart = sq(start + sw->loffset);
        vstop = sq(start + len - sw->roffset - 1);
        if (vstop < vstart)
         { fprintf(stderr,"Searching from %ld to %ld\n",vstart,psmax);
           fprintf(stderr,"Searching from 1 to %ld\n",vstop); }
        else
         fprintf(stderr,"Searching from %ld to %ld\n",vstart,vstop); }
     if (sw->both != 1)
      { sw->start = start;
        sw->comp = 0;
        nt = find_mt_trna(d,seq,len,nt,sw); }
     if (sw->both > 0)
      { sense_switch(seq,cseq,len);
        sw->start = start+len;
        sw->comp = 1;
        nt = find_mt_trna(d,cseq,len,nt,sw); }
     if (!flag)
      { s = seq;
        sf = se - drewind;
        se = seq + drewind;
        while (s < se) *s++ = *sf++;
        start += len - drewind;
        goto NX; }
     disp_gene_set(d,ts,nt,sw);
     if (sw->verbose) fprintf(stderr,"%s\nSearch Finished\n\n",d->seqname);
     d->ns++; }
  if (d->ns > 1)
   { fprintf(f,"\n\n%d sequences searched\n",d->ns);
     if (sw->trna | sw->mtrna)
      { fprintf(f,"Total tRNA genes = %d\n",sw->ngene[0]);
	    if (sw->matchacceptor)
         fprintf(f,"Total iso-acceptor mismatches = %d\n",sw->iamismatch); }}
  free((void *)ts); }


void bopt_fastafile(data_set *d, csw *sw)
{ int i,nt,flag,len;
  int *s,*sf,*se,*sc,*swrap;
  int seq[2*LSEQ+WRAP+1],cseq[2*LSEQ+WRAP+1],wseq[2*WRAP+1];
  long gap,start,rewind,drewind,psmax,tmaxlen,vstart,vstop;
  FILE *f = sw->f;
  ts = (gene *)malloc(NT*sizeof(gene));
  if (ts == NULL)
   { fprintf(stderr,"Not enough memory available to store detected genes\n");
     exit(1); }
  sw->space = NT;
  rewind = MAXTAGDIST + 20;
  if (sw->trna | sw->mtrna)
   { tmaxlen = MAXTRNALEN + sw->maxintronlen;
     if (rewind < tmaxlen) rewind = tmaxlen; }
  sw->loffset = rewind;
  sw->roffset = rewind;
  drewind = 2*rewind;
  d->ns = 0;
  d->nextseq = 0L;
  while (d->nextseq >= 0L)
   { d->seqstart = d->nextseq;
     if (!seq_init(d)) break;
     psmax = d->psmax;
     if (sw->verbose)
      { fprintf(stderr,"%s\n",d->seqname);
        fprintf(stderr,"%ld nucleotides in sequence\n",psmax);
        fprintf(stderr,"Mean G+C content = %2.1f%%\n",100.0*d->gc); }
     if (sw->batch < 2) fprintf(f,">%s\n",d->seqname);
     init_gene(ts,0,NT);
     nt = 0;
     flag = 0;
     start = 1L;
     se = seq;
     if (sw->linear)
      { for (i = 0; i < rewind; i++) *se++ = NOBASE;
        start -= rewind; }
     else
      { if (psmax <= drewind)
         { gap = drewind - psmax;
           sc = se + gap;
           while (se < sc) *se++ = NOBASE;
           swrap = wseq;
           sc = se + psmax;
           while (se < sc)
            { *se = move_forward(d);
              *swrap++ = *se++; }
           sc = swrap + gap;
           while (swrap < sc) *swrap++ = NOBASE;
                   swrap = wseq;
                   sc = swrap + psmax;
                   while (swrap < sc) *se++ = *swrap++;
                   swrap = wseq;
                   sc = swrap + drewind;
                   while (swrap < sc) *se++ = *swrap++;
                   sw->loffset = drewind;
                   sw->roffset = drewind;
                   start -= drewind;
                   flag = 1;
                   goto SH; }
        else
         { swrap = wseq;
           sc = seq + drewind;
           while (se < sc)
            { *se = move_forward(d);
              *swrap++ = *se++; }}}
     sc = seq + LSEQ;
     NX:
     while (se < sc)
      { *se++ = move_forward(d);
        if (d->ps >= psmax)
         { if (sw->linear)
            for (i = 0; i < rewind; i++) *se++ = NOBASE;
           else
            { sc = wseq + drewind;
              swrap = wseq;
              while (swrap < sc) *se++ = *swrap++; }
           flag = 1;
           break; }}
         SH:
     len = (int)(se - seq);
     if (sw->verbose)
      { vstart = sq(start + sw->loffset);
        vstop = sq(start + len - sw->roffset - 1);
        if (vstop < vstart)
         { fprintf(stderr,"Searching from %ld to %ld\n",vstart,psmax);
           fprintf(stderr,"Searching from 1 to %ld\n",vstop); }
        else
         fprintf(stderr,"Searching from %ld to %ld\n",vstart,vstop); }
     if (sw->both != 1)
      { sw->start = start;
        sw->comp = 0;
        nt = find_mt_trna(d,seq,len,nt,sw); }
     if (sw->both > 0)
      { sense_switch(seq,cseq,len);
        sw->start = start+len;
        sw->comp = 1;
        nt = find_mt_trna(d,cseq,len,nt,sw); }
     if (!flag)
      { s = seq;
     sf = se - drewind;
     se = seq + drewind;
     while (s < se) *s++ = *sf++;
     start += len - drewind;
     goto NX; }
     batch_gene_set(d,ts,nt,sw);
     if (sw->verbose) fprintf(stderr,"%s\nSearch Finished\n\n",d->seqname);
     d->ns++; }
  if ((d->ns > 1) && (sw->batch < 2))
   { fprintf(f,">end \t%d sequences",d->ns);
     if (sw->trna || sw->mtrna) fprintf(f," %d tRNA genes",sw->ngene[0]);
     fputc('\n',f); }
  free((void *)ts); }


void arwen_help_menu()
{ printf("\n");
  printf("---------------------------\n");
  printf("ARWEN v1.2 Dean Laslett\n");
  printf("---------------------------\n");
  printf("\n");
  printf("Please reference the following paper if you use this\n");
  printf("program as part of any published research.\n\n");
  printf("Laslett, D. and Canback, B. (2008) ARWEN: a\n");
  printf("program to detect tRNA genes in metazoan mitochondrial\n");
  printf("nucleotide sequences\n");
  printf("Bioinformatics, 24(2); 172-175.\n\n\n");
  printf("ARWEN detects metazoan mitochondrial tRNA genes.\n");
  printf("\n");
  printf("Usage:\n");
  printf("arwen -v -s -d -c -l -a -w -j -seq -mt");
  printf(" -tv -gc -o <outfile> <filename>\n\n");
  printf("<filename> is assumed to contain one or more sequences\n");
  printf("in FASTA format. Results of the search are printed to\n");
  printf("STDOUT. All switches are optional and case-insensitive.\n");
  printf("\n");
  printf("    -mtmam        Search for Mammalian mitochondrial tRNA\n");
  printf("                  genes. -tv switch set.\n");
  printf("                  Mammalian mitochondrial genetic code used.\n");
  printf("    -mtx          Low scoring tRNA genes are\n"); 
  printf("                  not reported.\n");
  printf("    -gc<num>      Use the GenBank transl_table = <num> genetic code.\n");
  printf("    -gcstd        Use standard genetic code.\n");
  printf("    -gcmet        Use composite Metazoan mitochondrial genetic code.\n");
  printf("    -gcvert       Use Vertebrate mitochondrial genetic code.\n");
  printf("    -gcinvert     Use Invertebrate mitochondrial genetic code.\n");
  printf("    -gcyeast      Use Yeast mitochondrial genetic code.\n");
  printf("    -gcprot       Use Mold/Protozoan/Coelenterate");
  printf(" mitochondrial genetic code.\n");
  printf("    -gcciliate    Use Ciliate genetic code.\n");
  printf("    -gcflatworm   Use Echinoderm/Flatworm mitochondrial genetic code.\n");
  printf("    -gceuplot     Use Euplotid genetic code.\n");
  printf("    -gcbact       Use Bacterial/Plant Chloroplast genetic code.\n");
  printf("    -gcaltyeast   Use alternative Yeast genetic code.\n");
  printf("    -gcascid      Use Ascidian Mitochondrial genetic code.\n");
  printf("    -gcaltflat    Use alternative Flatworm Mitochondrial genetic code.\n");
  printf("    -gcblep       Use Blepharisma genetic code.\n");
  printf("    -gcchloroph   Use Chlorophycean Mitochondrial genetic code.\n");
  printf("    -gctrem       Use Trematode Mitochondrial genetic code.\n");
  printf("    -gcscen       Use Scenedesmus obliquus Mitochondrial genetic code.\n");
  printf("    -gcthraust    Use Thraustochytrium Mitochondrial genetic code.\n");
  printf("                  Individual modifications can be appended using\n");
  printf("    ,BBB=<aa>     B = A,C,G, or T. <aa> is the three letter\n");
  printf("                  code for an amino-acid. More than one modification\n");
  printf("                  can be specified. eg -gcvert,aga=Trp,agg=Trp uses\n");
  printf("                  the Vertebrate Mitochondrial code and the codons\n");
  printf("                  AGA and AGG changed to Tryptophan.\n");            
  printf("    -tv           Do not search for mitochondrial ");
  printf("TV replacement\n");
  printf("                  loop tRNA genes.\n");
  printf("    -c            Assume that each sequence has a circular\n");
  printf("                  topology. Search wraps around each end.\n");
  printf("                  Default setting.\n");
  printf("    -l            Assume that each sequence has a linear\n");
  printf("                  topology. Search does not wrap.\n");
  printf("    -d            Double. Search both strands of each\n");
  printf("                  sequence. Default setting.\n");
  printf("    -s or -s+     Single. Do not search the complementary\n");
  printf("                  (antisense) strand of each sequence.\n");
  printf("    -sc or -s-    Single complementary. Do not search the sense\n");
  printf("                  strand of each sequence.\n");
  printf("    -rp           Report low scoring genes as possible pseudogenes\n");
  printf("    -seq          Print out primary sequence.\n");
  printf("    -br           Show secondary structure of tRNA gene primary\n");
  printf("                  sequence with round brackets.\n");
  printf("    -fasta        Print out primary sequence in fasta format.\n");
  printf("    -fo           Print out primary sequence in fasta format only\n");
  printf("                  (no secondary structure).\n");
  printf("    -j            Display 4-base sequence on 3' end of astem\n");
  printf("                  regardless of predicted amino-acyl acceptor\n");
  printf("                  length.\n");
  printf("    -jr           Allow some divergence of 3' ");
  printf("amino-acyl acceptor\n");
  printf("                  sequence from NCCA.\n");
  printf("    -jr4          Allow some divergence of 3' ");
  printf("amino-acyl acceptor\n");
  printf("                  sequence from NCCA, and display 4 bases.\n");
  printf("    -v            Verbose. Prints out search progress\n");
  printf("                  to STDERR.\n");
  printf("    -o <outfile>  print output into <outfile>. If <outfile>\n");
  printf("                  exists, it is overwritten.\n");
  printf("                  By default, output goes to STDOUT.\n");
  printf("    -w            Print out genes in batch mode.\n");
  printf("                  Output is in the form:\n\n");
  printf("                  Sequence name\n");
  printf("                  N genes found\n");
  printf("                  1 tRNA-<species> [locus 1]");
  printf(" <Apos> (nnn)\n");
  printf("                            .          \n");
  printf("                            .          \n");
  printf("                  N tRNA-<species> [Locus N]");
  printf(" <Apos> (nnn)\n");
  printf("\n                  N is the number of genes found\n");
  printf("                  <species> is the tRNA iso-acceptor species\n");
  printf("                  <Apos> is the tRNA anticodon ");
  printf("relative position\n");
  printf("                  (nnn) is the tRNA anticodon base triplet\n");
  printf("\n\n"); }

void error_report(int n, char *s)
{ switch(n)
   { case 0: fprintf(stderr,
              "-%s not recognised, type arwen -h for help\n",s);
             break;
     case 1: fprintf(stderr,
              "-%s not understood, type arwen -h for help\n",s);
             break;
     case 2: fprintf(stderr,"Could not open %s\n",s);
             break;
     case 3: fprintf(stderr,
             "No sequence file specified, type arwen -h for help\n");
             break;
     case 4: fprintf(stderr,"Don't know genetic code %s\n",s);
             break;
     case 5: fprintf(stderr,"Too many genetic code modifications (max=%d)\n",
                     MAXGCMOD);
             break;
     default: break; }
  exit(0); }


void process_genecode_switch(char *s, csw *sw)
{ int i,m,lmax,len[NGENECODE],anticodon,b[3];
  char c,*ss,*se;
  static char genecodetag[NGENECODE][10] =
   { "MET",
     "STD","VERT","YEAST","PROT","INVERT",
     "CILIATE","DELETED","DELETED","FLATWORM","EUPLOT",
     "BACT","ALTYEAST","ASCID","ALTFLAT","BLEP",
     "CHLOROPH","DELETED","DELETED","DELETED","DELETED",
     "TREM","SCEN","THRAUST" };
  sw->geneticcode = STANDARD;
  sw->gcfix = 1;
  c = *s;
  if (c >= '0')
   if (c <= '9')
    { iconvert(s,&i);
      if ((i >= 0) && (i < NGENECODE)) sw->geneticcode = i;
      goto MOD; }
  for (i = 0; i < NGENECODE; i++) 
   { len[i] = 0;
     ss = s;
     se = genecodetag[i];
     while (c == *ss++)
      { if (upcasec(c) != *se++) break;
        len[i]++; }}
  m = -1;
  lmax = 0;
  i = -1;
  while (++i < NGENECODE)
   if (len[i] > lmax)
    { m = i;
      lmax = len[i]; }
  if (m >= 0) sw->geneticcode = m;
  else error_report(4,s);
  MOD:
  sw->ngcmod = 0;
  ss = s;
  while (ss = strpos(ss,","))
   { if (sw->ngcmod >= MAXGCMOD) error_report(5,NULL);
     ss++;
     for (i = 0; i < 3; i++)
      { b[i] = Adenine;
        c = upcasec(ss[i]);
        if (c == 'C') b[i] = Cytosine;
        if (c == 'G') b[i] = Guanine;
        if (c == 'T') b[i] = Thymine;
        if (c == 'U') b[i] = Thymine; }
      anticodon = ((Thymine - b[2])<<4) + ((Thymine - b[1])<<2) + (Thymine - b[0]);
     if (!(se = strpos(ss,"="))) break;
     se++;
     for (i = 0; i < NAMINOACID; i++)
      if (upcasec(se[0]) == upcasec(aaname[i][0]))
       if (upcasec(se[1]) == upcasec(aaname[i][1]))
        if (upcasec(se[2]) == upcasec(aaname[i][2]))
         { aamap[sw->geneticcode][anticodon] = i;
           sw->gcmod[sw->ngcmod] = anticodon;
           sw->ngcmod++;
           break; } }}
    

int main(int z, char *v[])
{ int i,l,lv,filecounter;
  char c1,c2,c3,*s;
  double psthresh;
  data_set d;
  static csw sw =
   { "",NULL,0,0,1,1,1,0,
     METAZOAN_MT,0,{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},0,METAZOAN_MT,
     1,0,5,5,1,0,0,0,2,0,0,0,0,0,0,2,1,0,
     132.0,325.0,10.85,4.0,29.0,26.0,7.5,14.0,10.0,25.0,
     1,9,8,0,0,0,0,mtRNAtthresh,mtRNAdthresh,mtRNAdtthresh,-7.9,-6.0,
     {0,0,0,0,0},0,0,0,0,0,0,
     { 45, 45, 45, 45, 45, 45, 45, 45, 45, 45,
       45, 45, 45, 45, 45, 45, 45, 45, 45, 45,
       45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 10,
       65, 82, 87, 69, 78, 32,
       118, 49, 46, 50, 32, 32, 32, 68, 101, 97,
       110, 32, 76, 97, 115, 108, 101, 116, 116, 10,
       45, 45, 45, 45, 45, 45, 45, 45, 45, 45,
       45, 45, 45, 45, 45, 45, 45, 45, 45, 45,
       45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 10, TERM }};
  sw.f = stdout;
  filecounter = 0;
  i = 0;
  while (++i < z)
   if (*(v[i]) == '-')
    { lv = length(v[i]);
      if (lv < 2) continue;
      s = v[i] + 1;
      c1 = upcasec(*s);
      c2 = (lv > 2)?upcasec(s[1]):' ';
      c3 = (lv > 3)?upcasec(s[2]):' ';
      switch(c1)
       { case  'E': sw.energydisp = (c2 == 'S')?2:1;
                    break;
         case  'A': if (c2 == '7') sw.extastem = 0;
                    else
                     if (c2 == 'A') sw.matchacceptor = 1;
                    break;
         case  'B': if (c2 == 'R') sw.seqdisp = 2;
                    else sw.libflag = 1;
                    break;
         case  'X': sw.libflag = 2;
                    break;
         case  'W': if (sw.batch < 1) sw.batch = 1;
                    break;
         case  'V': sw.verbose = 1;
                    break;
         case  'S': if (c2 == 'E')
                     { if (sw.seqdisp < 1) sw.seqdisp = 1;
                       break; }
                    if ((c2 == 'C') || (c2 == '-'))
                     { sw.both = 1;
                       break; }
                    sw.both = 0;
                    break;
         case  'F': if (softstrpos(s,"O")) 
                     { sw.batch = 2;
                       if (softstrpos(s,"S")) sw.batch |= 0x4;
                       if (softstrpos(s,"N")) sw.batch |= 0x8;
                       if (softstrpos(s,"C")) sw.batch |= 0x10; }
                    else 
                     { if (softstrpos(s,"C")) sw.seqdisp = 4;
                       else sw.seqdisp = 3; }
                    break;
         case  'D': sw.both = 2;
                    break;
         case  'L': sw.linear = 1;
                    break;
         case  'C': if (c2 == '7') sw.cloop7 = 1;
                    else sw.linear = 0;
                    break;
         case  'J': if (lv > 2)
                     { if (c2 == 'R') sw.aataildiv = 1;
                       if (c3 == '4') sw.aataildisp = 1; }
                    else sw.aataildisp = 1;
                    break;
         case  'T': if (c2 == 'V')
                     { sw.tvloop = 0;
                       break; }
                    break;
         case  'M': if (c2 == 'T')
                     { if (!sw.gcfix) sw.geneticcode = METAZOAN_MT;
                       if (lv > 3)
                        { s += 2;
                          c3 = upcasec(*s);
                          if (c3 == 'M')
                           { do c3 = upcasec(*++s);
                             while ((c3 == 'A') || (c3 == 'M')
                                     || (c3 == 'L'));
                             sw.tvloop = 0;
                             sw.geneticcode = VERTEBRATE_MT;
                             sw.discrim = MAMMAL_MT; }
                          if (c3 == 'X')
                           { c3 = upcasec(*++s);
                             sw.mtxdetect = 0; }
                          if (c3 != '-')
                           if (c3 != '.')
                            if ((c3 < '0') || (c3 > '9'))
                             break;
                          s = dconvert(s,&sw.mtdtthresh);
                          if (*s == ',') s = dconvert(s+1,&sw.mttthresh);
                          if (*s == ',') s = dconvert(s+1,&sw.mtdthresh);
                          if (*s == ',') s = dconvert(s+1,&sw.mttarmthresh);
                          if (*s == ',') dconvert(s+1,&sw.mtdarmthresh); }}
                    break;
         case  'P': if (c2 == 'S')
                     { if (c3 != '-')
                        if (c3 != '.')
                         if ((c3 < '0') || (c3 > '9'))
                          { sw.mtdtthresh *= PSEUDOGENElevel;
                            sw.mttthresh *= PSEUDOGENElevel;
                            sw.mtdthresh *= PSEUDOGENElevel;
                            break; }
                       psthresh = 1.0;
                       dconvert(s+2,&psthresh);
                       sw.mtdtthresh *= psthresh;
                       sw.mttthresh *= psthresh;
                       sw.mtdthresh *= psthresh;
                       break; }
                    break;
         case  'G': if (c2 != 'C') break;
                    process_genecode_switch(s+2,&sw);
                    break;
         case  'R': if (c2 == 'N') sw.repeatsn = 1;
                    else 
                     if (c2 == 'P') sw.reportpseudogenes = 1;
                    break;
         case  'Q': sw.showconfig = 0;
                    break;
         case  'H': arwen_help_menu();
                    exit(0);
         case  'O': if (lv > 2)
                     s++;
                    else
                     { if (++i >= z) break;
                       s = v[i]; }
                    sw.f = fopen(s,"w");
                    if (!sw.f) error_report(2,s);
                    break;
         default:   error_report(0,s); }}
   else
    if (filecounter < 1)
     { d.f = fopen(v[i],"r");
       if (d.f)
        filecounter++;
       else
        error_report(2,v[i]); }
    else
     if (filecounter < 2)
      { sw.f = fopen(v[i],"w");
        if (!sw.f) error_report(2,v[i]);
        filecounter++; }
     else
     error_report(0,v[i]);
  if (filecounter < 1)
   error_report(3,NULL);
  if (sw.mtrna) sw.trna = 0;
  if (sw.libflag) fprintf(sw.f,"tRNA Library\n");
  if (sw.batch) bopt_fastafile(&d,&sw);
  else iopt_fastafile(&d,&sw);
  fclose(d.f);
  if (!sw.batch && sw.showconfig)
   { fprintf(sw.f,"Configuration: ");
     i = -1;
     while (++i < z) fprintf(sw.f,"%s ",v[i]);
     fputc('\n',sw.f); }
  if (sw.f != stdout) fclose(sw.f);
  return(0); }

