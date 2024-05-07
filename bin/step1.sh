export PATH=/home/kino/kino/kit/atat/bin:$PATH
#extract_vasp
#cp str_relax.out str.out
python /home/kino/kino/kit/python3atat/bin/poscar2str.py opt.vasp
cp str.out str_relax.out
fitfc -er=12 -ns=1 -dr=0.1 -nrr

for name in vol_0/p*; 
do (echo $name; cd $name;
# python /home/kino/kino/kit/python3atat/bin/str2poscar.py str.out 
/home/kino/kino/kit/python3atat/bin/str2poscar.sh
python /home/kino/kino/kit/python3atat/bin/atatposcarfix.py POSCAR
);
done

