# Unix
## Example of basic commands

list all contents in the current directory
```
ls
```
list all contents with details in the current directory
```
ls -l
```

list all contents that end with .txt
```
ls *.txt
```

create a new directory called "newDirectory"
```
mkdir newDirectory
```

change the current working directory to "newDirectory"
```
cd newDirectory
```

list contents of the directory above the current directory without being in that directory
```
ls ..
```

create a new file
```
> newfile01.txt
```
 or

 ```
touch newfile02.txt
```

or

```
echo "TEXT03" > newfile03.txt
```

remove a file
```
rm newfile01.txt
```

move a file
```
mv file02.txt newDirectory
```


copy a file
```
cp newfile01.txt newDirectory
```


navigate up one level in the directory structure
```
cd ..
```

remove directory (if it is empty)
```
rmdir newDirectory
```

print working directory
```
pwd
```

copy file "newfile.txt" to "copyfile.txt"
```
cp newfile.txt copyfile.txt
```

move file "newfile.txt" from current directory to "Dir" directory
```
mv newfile.txt Dir/
```

remove a file "copyfile.txt"
```
rm copyfile.txt
```

display contents of a file
```
cat newfile.txt 
```

print the first top 10 lines of a file
```
head newfile.txt
```

show data at the end of a file
```
tail newfile.txt
```

text search the word "pattern" in a file
```
grep "pattern" newfile02.txt
```

download and save a file from www to working directory
```
wget https://bioinformaticsworkbook.org/Appendix/Unix/assets/saltandpepper.jpg
```


## Example of installing and updating software

update ubuntu packages
```
sudo apt-get update
```

upgrade packages
```
sudo apt-get upgrade
```

install a package in ubuntu
```
sudo apt-get install packageName
```


## References
https://www.doc.ic.ac.uk/~wjk/UnixIntro/Lecture3.html

https://users.cs.duke.edu/~alvy/courses/unixtut/unix1.html

https://datacarpentry.org/2015-11-04-ACUNS/shell-intro/

https://bioinformaticsworkbook.org/Appendix/Unix/unix-basics-1.html#gsc.tab=0
