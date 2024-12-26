# Basic Linux Commands for Biology Students

## Table of Contents
1. [Introduction](#introduction)  
2. [Navigating Directories](#navigating-directories)  
3. [Viewing & Manipulating Files](#viewing--manipulating-files)  
4. [Searching Files & Directories](#searching-files--directories)  
5. [User Permissions](#user-permissions)  
6. [Working with Processes](#working-with-processes)  
7. [Wrap-Up](#wrap-up)

---
<br>

## Introduction

In biology research, Linux often underpins servers and computing clusters used for tasks like sequencing analysis, protein structure modeling, or data-intensive pipelines. Understanding basic commands saves time, reduces errors, and lets you automate repetitive work.

### Opening a Terminal
- **Linux machine**: You can open a terminal (shell) directly.  
- **Windows**: You might use [Windows Subsystem for Linux (WSL)](https://learn.microsoft.com/en-us/windows/wsl/install) or tools like PuTTY for SSH.  
- **macOS**: Use the built-in **Terminal** application.  

### Goals
- Familiarize yourself with the core commands.  
- Practice applying these commands to real-world biological data scenarios.  
- Build confidence in navigating and manipulating files in the Linux environment.

---

### Quiz 1 (True/False?)
1. **T/F**: Linux is rarely used in research computing environments.  
2. **T/F**: You can open a Linux terminal in macOS without installing anything extra.  
3. **T/F**: Knowing basic Linux commands can help automate repetitive tasks in biology research.

<details>
  <summary>Click to Reveal Answers</summary>

1. **False** – Linux is *very* commonly used in research.  
2. **True** – macOS has a built-in Terminal that supports many UNIX commands.  
3. **True** – Automation is a major benefit of command-line usage.
</details>

---
<br>

## Navigating Directories

### Core Concepts
1. **pwd** – Stands for “print working directory” (shows where you are).  
2. **ls** – Lists the contents of a directory. Common options:  
   - **ls -l** (long listing, shows file permissions and sizes)  
   - **ls -a** (shows hidden files)  
3. **cd** – Changes the current directory. Examples:  
   - **cd /home/student/data** (go to `/home/student/data`)  
   - **cd ..** (move up one directory level)  
   - **cd ~** (go to your home directory)  

### Exercise
1. Open your terminal.  
2. Use **pwd** to see your current directory.  
3. Use **ls** to view its contents.  
4. Move one directory up with **cd ..**, then list contents again.

### Quiz 2 (True/False?)
1. **T/F**: `cd ..` takes you to the parent (upper) directory.  
2. **T/F**: `pwd` shows you the contents of your current directory.  
3. **T/F**: `ls -a` displays hidden files as well.

<details>
  <summary>Click to Reveal Answers</summary>

1. **True**  
2. **False** – `pwd` shows the path of the current directory, not the contents.  
3. **True**
</details>

---
<br>

## Viewing & Manipulating Files

### Core Concepts
1. **cat** – Displays the contents of a file in the terminal.  
2. **head / tail** – Shows the first or last N lines of a file (default is 10 lines). Examples:  
   - `head sample.txt` (first 10 lines)  
   - `tail -n 5 sample.txt` (last 5 lines)  
3. **mkdir** – Creates a new directory.  
4. **rm** – Removes files or directories (use with caution!).  
   - `rm file.txt` deletes `file.txt`.  
   - `rm -r folder_name` removes a folder and all its contents recursively.  
5. **cp** – Copies a file or directory.  
   - `cp old.txt new.txt` creates a copy.  
   - `cp -r old_folder new_folder` copies folders.  
6. **mv** – Moves or renames a file or directory.  
   - `mv file.txt /home/student/data/` (move file)  
   - `mv old.txt new.txt` (rename file)  

### Exercise
1. Create a new directory called **test_dir** (`mkdir test_dir`).  
2. Create a small text file with some lines (e.g., `nano sample.txt`, add text, save).  
3. Use **head** or **tail** to view part of the file.  
4. Copy `sample.txt` into **test_dir** using `cp sample.txt test_dir/`.  
5. Remove `sample.txt` with `rm sample.txt`.

### Quiz 3 (True/False?)
1. **T/F**: `head -n 3 sample.txt` will show the last 3 lines of `sample.txt`.  
2. **T/F**: `cp -r` is used to copy directories.  
3. **T/F**: `rm` can delete directories without any special flags.

<details>
  <summary>Click to Reveal Answers</summary>

1. **False** – `head -n 3` shows the *first* 3 lines.  
2. **True**  
3. **False** – You need `-r` (recursive) to delete directories with `rm`.
</details>

---
<br>

## Searching Files & Directories

### Core Concepts
1. **grep** – Searches for text patterns in files. Examples:  
   - `grep "ATG" sequence.fasta` (search for the string “ATG”)  
   - `grep -i "gene" results.txt` (case-insensitive search for “gene”)  
2. **find** – Locates files and directories based on criteria like name or size. Examples:  
   - `find . -name "*.fasta"` (search recursively from current directory for files ending with `.fasta`)  
   - `find /home/student/data -type d` (lists all directories in `/home/student/data`)  

### Exercise
1. Create a file containing some text lines (e.g., `nano biology_terms.txt`) with words like “cell”, “Gene”, “GENE”, “DNA”, etc.  
2. Use `grep -i "gene" biology_terms.txt` and observe how many lines contain “gene” in any combination of uppercase/lowercase.  
3. Use `find . -name "*terms.txt"` to locate the file you just created.

### Quiz 4 (True/False?)
1. **T/F**: `grep -i "DNA" file.txt` will search for “DNA” only in uppercase, ignoring lower or mixed case.  
2. **T/F**: `find . -name "*.txt"` will look for files in the current directory (and subdirectories) ending with `.txt`.  
3. **T/F**: Using `grep "genome"` without specifying a file will print an error message.

<details>
  <summary>Click to Reveal Answers</summary>

1. **False** – The `-i` option does a case-insensitive search (so it will match “dna”, “DNA”, “Dna”, etc.).  
2. **True**  
3. **True** – You must specify a file or a stream for `grep`.
</details>

---
<br>

## User Permissions

### Core Concepts
1. **ls -l** Output Format:  
   ```
   drwxr-xr-x  2 user group 4096 Aug 30 12:34 folder
   -rw-r--r--  1 user group  256 Aug 30 12:34 file.txt
   ```
   - **d** or `-` at the start indicates directory (d) or file (-).  
   - Next 9 characters are permissions for **owner**, **group**, and **others** (e.g., `rwxr-xr-x`).  
2. **chmod** – Changes a file or directory's permissions. Examples:  
   - `chmod 755 script.sh` (owner can read/write/execute, group & others can read/execute).  
3. **chown** – Changes a file's owner or group (requires elevated privileges).

### Exercise
1. Run `ls -l` on some directory. Note the permission bits (e.g., `rw-r--r--`).  
2. Try changing permissions for a file you own (`chmod 600 file.txt`) and see how they change in `ls -l`.  

### Quiz 5 (True/False?)
1. **T/F**: `drwxr--r--` indicates a directory in which the owner can read, write, execute, but group and others can only read.  
2. **T/F**: `chmod` is used to change file ownership.  
3. **T/F**: `ls -l` displays user permissions for each file/directory.

<details>
  <summary>Click to Reveal Answers</summary>

1. **True** – `d` means directory, `rwx` for owner, `r--` for group, `r--` for others.  
2. **False** – `chmod` changes permissions, while `chown` changes ownership.  
3. **True**  
</details>

---
<br>

## Working with Processes

### Core Concepts
1. **ps** – Lists current running processes.  
2. **top** / **htop** – Provides a dynamic, real-time view of running processes, CPU usage, and memory usage.  
3. **kill** – Sends a signal to a process (commonly to terminate it). Example:
   - `kill 12345` will terminate the process with the PID (Process ID) 12345.
4. **&** – Run a command in the background. Example:
   - `python long_script.py &`  

### Exercise
1. In one terminal, run something like `sleep 100 &` which simulates a process running in the background for 100 seconds.  
2. In another terminal, run `ps` or `top` to see that process.  
3. Use `kill` to end that process. Check with `ps` again to confirm it ended.

### Quiz 6 (True/False?)
1. **T/F**: `ps` shows a static list of running processes.  
2. **T/F**: `kill 0` kills all processes on the system.  
3. **T/F**: Adding `&` to the end of a command runs it in the background.

<details>
  <summary>Click to Reveal Answers</summary>

1. **True** – `ps` outputs a one-time snapshot. (Whereas `top` is dynamic.)  
2. **False** – `kill 0` won’t kill *all* processes (that would require more advanced usage or root privileges).  
3. **True**
</details>

---
<br>

## Wrap-Up

You have learned how to:
1. Navigate directories (`pwd`, `ls`, `cd`)  
2. Manipulate files (`mkdir`, `rm`, `cp`, `mv`)  
3. Search files (`grep`, `find`)  
4. Understand file permissions (`ls -l`, `chmod`)  
5. Manage processes (`ps`, `kill`, `top`)  

### Additional Resources
- **Linux Command Line Tutorial** by [Linux Journey](https://linuxjourney.com/)  
- **Software Carpentry** lessons: [Unix Shell](https://software-carpentry.org/lessons/)  
