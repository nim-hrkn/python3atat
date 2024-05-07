export PATH=/home/kino/kino/kit/atat/bin:$PATH
cp energy vol_0
cp str_relax.out vol_0

for name in vol_0/p*; 
do (echo $name; cd $name; 
python /home/kino/kino/kit/python3atat/bin/poscar2str.py --strout str_relax.out POSCAR; 
python /home/kino/kino/kit/python3atat/bin/extract_force.py opt.json );
done

# fitfc -f -frnn=10  -fn -dT=10
fitfc -f -fr=10  -fn -dT=10
# fitfc -f -frnn=10 -df=kpath -fn
fitfc -f -fr=10 -df=kpath -fn

