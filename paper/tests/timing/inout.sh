#/bin/bash
Cx=.2
Cy=.5
nt=64
it=3

#for it in {1,2,3}; do
  for nxy in {0064,0096,0128,0192,0256,0384,0512,0768,1024,1536,2048,3072,4096}; do
    ../../cpp/inout-cpp $nxy $nxy $Cx $Cy $nt $it \
      1>"nx=${nxy}_ny=${nxy}_Cx=${Cx}_Cy=${Cy}_nt=${nt}_it=${it}:in" \
      2>"nx=${nxy}_ny=${nxy}_Cx=${Cx}_Cy=${Cy}_nt=${nt}_it=${it}:out"
  done
#done
