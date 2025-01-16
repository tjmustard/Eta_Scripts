for file in *.txt
do
    base=`basename $file .txt`
    echo "#CONFIG" > $base.tm
    cat $file >> $base.tm
    echo "#END" >> $base.tm
done
