#/bin/bash
Cx=.2
Cy=.5
it=3
nt=100

for nxy in {0032,0064,0128,0256,0512,1024}; do
  ../../cpp/inout-cpp $nxy $nxy $Cx $Cy $nt $it \
    1>"nx=${nxy}_ny=${nxy}_Cx=${Cx}_Cy=${Cy}_nt=${nt}_it=${it}:in" \
    2>"nx=${nxy}_ny=${nxy}_Cx=${Cx}_Cy=${Cy}_nt=${nt}_it=${it}:out"
done
