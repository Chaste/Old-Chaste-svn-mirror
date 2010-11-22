# Created by Chaste version 2.1.10890 on Mon, 22 Nov 2010 03:19:59 +0000.  Chaste was built on Mon, 22 Nov 2010 03:19:47 +0000 by machine (uname) 'Linux compphys11 2.6.32-25-generic #45-Ubuntu SMP Sat Oct 16 19:52:42 UTC 2010 x86_64' using settings: default, no Chaste libraries.
gfx read node SimulationResults.exnode 
gfx read elem SimulationResults.exelem generate_faces_and_lines 
 gfx cre win 1 
for ($i=0; $i<2; $i++) { 
  gfx read node cube_2mm_12_elements_$i.exnode time $i
}
