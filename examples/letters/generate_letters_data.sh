if [ ! -e ../../build/strokestrip ]
then
  echo "../../build/strokestrip does not exist. Have you built the code yet? (See the Building section of README.md)"
  exit 1
fi

for f in *.scap
do
  echo "Running on $f"
  ../../build/strokestrip "$f"
done
