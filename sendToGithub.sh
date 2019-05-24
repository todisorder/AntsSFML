echo "Inserir comentários de git commit"
read COMM
git add -A
git commit -m "$COMM"
git push

