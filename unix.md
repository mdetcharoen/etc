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

create a new directory call "newDirectory"
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
grep "pattern" newfile.txt
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
