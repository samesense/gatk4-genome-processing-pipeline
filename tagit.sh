git tag --delete 1.0.0
git push --delete origin 1.0.0

git add -A
git commit -m 'Resource'
git tag -a 1.0.0 `git log -1 | head -1 | cut -f 2 -d ' '` -m "Resource"

git push origin chop
git push origin 1.0.0
