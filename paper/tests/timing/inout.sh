#/bin/bash
Cx=.2
Cy=.5
it=3

#for it in {1,2,3}; do
  for nxy in {0064,0081,0102,0128,0161,0203,0256,0323,0406,0512,0645,0813,1024,1290,1626,2048,2580,3251,4096}; do
    nt=`echo "2^16/$nxy"|bc`
echo $nt
    ../../cpp/inout-cpp $nxy $nxy $Cx $Cy $nt $it \
      1>"nx=${nxy}_ny=${nxy}_Cx=${Cx}_Cy=${Cy}_nt=${nt}_it=${it}:in" \
      2>"nx=${nxy}_ny=${nxy}_Cx=${Cx}_Cy=${Cy}_nt=${nt}_it=${it}:out"
  done
#done
