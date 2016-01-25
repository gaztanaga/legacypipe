list=$1

wget --save-cookies=ptf.txt -O /dev/null "http://irsa.ipac.caltech.edu/account/signon/login.do?josso_cmd=login&josso_username=kburleigh@lbl.gov&josso_password=Jasper88"
for i in `cat ${list}`
do
    wget --load-cookies=ptf.txt http://irsa.ipac.caltech.edu/ibe/data/ptf/images/level1/${i}
done

