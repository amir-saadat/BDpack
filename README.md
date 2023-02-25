Brownain Dynamics Package (BDpack)
==================================

NOTE: the master branch has recently been renamed to main. Please follow these steps on your local branch to reflect this change on your local repo.

### Switch to the "master" branch:
$ git checkout master

### Rename it to "main":
$ git branch -m master main

### Get the latest commits (and branches!) from the remote:
$ git fetch

### Remove the existing tracking connection with "origin/master":
$ git branch --unset-upstream

### Create a new tracking connection with the new "origin/main" branch:
$ git branch -u origin/main
