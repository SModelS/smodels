**create a new repository:**
git ini <repository>

**checkout a repository:**
git clone /path/to/repository
using the remote server: git clone username@host:/path/to/repository


**workflow: add & commit**
git add <repository>
git add *
git commit -m "Commit message"

**pushing changes:**
git push origin master
Remote server: git remote add origin <server>

**branching:**
create a new branch: git checkout -b <branch>
switch back to master: git checkout master
delete the branch: git branch -d <branch>
a branch is not available to others unless you push the branch to your remote repository: git push origin <branch>


**update & merge**
