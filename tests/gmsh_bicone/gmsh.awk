# usage: gawk -f gmsh.awk  file.msh > surface_msh.inp
BEGIN {inod=0;iel=0;inn=0;its=0}
 /\$Nodes/ {inod=1}
 /\$EndNodes/ {inod=0}
 /\$Elements/ {iel=1}
 /\$EndElements/ {iel=0}
 inod==1&&NF==4 { inn++;xn[inn]=$2; yn[inn]=$3; zn[inn]=$4}
 iel==1&&$2==2 {its++;xts[its]=$6; yts[its]=$7; zts[its]=$8}
END{
 print inn
 i=1
 while (i<=inn) {
  print xn[i],yn[i],zn[i]
  i++}
 print its
 i=1
 while (i<=its) {
  print xts[i],yts[i],zts[i]
  i++}
}
