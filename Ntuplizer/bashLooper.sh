texts=(Saba Sabaa Sabaaa Sabaaaa)
nLoops=(3 1 1 2)

for value in {0..3}
do
    echo text: ${texts[$value]}
    echo running the file for the $((value+1)) time.
    python pythonSample.py ${texts[$value]} ${nLoops[$value]}
done

echo All done!
