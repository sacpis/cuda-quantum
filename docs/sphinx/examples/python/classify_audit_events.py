# ============================================================================ #
# Copyright (c) 2022 - 2024 NVIDIA Corporation & Affiliates.                   #
# All rights reserved.                                                         #
#                                                                              #
# This source code and the accompanying materials are made available under     #
# the terms of the Apache License 2.0 which accompanies this distribution.     #
# ============================================================================ #

import sys

def audit_hook(event, args):
    print(f"Audit event: {event}, Arguments: {args}")

sys.addaudithook(audit_hook)

# File Operations

# 'open' Event

with open("example.txt", "w") as file:
    file.write("Hello, World!")

# Classification: Unsafe: Opening sensitive system files or directories, especially in write mode.

# 'os.remove' Event

import os
os.remove("example.txt")

# Classification: Unsafe: Removing critical system files or directories.

# 'os.mkdir' Event

os.mkdir("new_directory")

# Classification: Unsafe: Creating directories in system or restricted areas.

# 'os.rmdir' Event

os.rmdir("test_directory")

# Classification: Unsafe: Removing critical system directories.

# 'os.rename' Event

os.rename("old_name.txt", "new_name.txt")

# Classification: Unsafe: Renaming critical system files or directories.

# 'os.chdir' Event

os.chdir("/tmp")

# Classification: Unsafe: Changing to directories containing sensitive information.

# 'os.chmod' Event

os.chmod("example.txt", 0x755)

# Classification: Unsafe: Changing permissions that could lead to security vulnerabilities.

# 'os.walk' Event

for root, dir, files in os.walk("."):
    for name in files:
        print(name)

# Classification: Unsafe: Walking through directories containing sensitive information.

# 'os.scandir' Event

for entry in os.scandir("."):
    print(entry.name)

# Classification: Unsafe: Scanning through directories containing sensitive information.

# 'shutil.copyfile' Event

import shutil

shutil.copyfile("source.txt", "destination.txt")

# Classification: Unsafe: Copying sensitive files to untrusted locations.

# 'os.utime' Event

os.utime("example.txt", None)

# Classification: Unsafe: Updating timestamps of sensitive or system files.

# Process Management

# 'os.system' Event

os.system("ls -la")

# Classification: Unsafe: Running commands that can modify system state, delete files, or execute arbitrary code

# 'subprocess.Popen' Event

import subprocess

subprocess.Popen(["echo", "Hello, World!"])

# Classification: Unsafe: Running subprocesses with user-supplied or untrusted commands.

# 'os.fork' Event

if os.fork() == 0:
    print("Child process")
else:
    print("Parent process")

# Classification: Unsafe: Forking processes that may lead to resource exhaustion or uncontrolled behavior.

# 'os.forkpty' Event

pid, fd = os.forkpty()
if pid == 0:
    print("Child process")
else:
    print("Parent process")

# Classification: Unsafe: Forking processes that may lead to resource exhaustion or uncontrolled behavior.

# 'os.exec' Event

os.execv("/bin/echo", ["echo", "Hello, World!"])

# Classification: Unsafe: Executing commands that can alter system state or execute arbitrary code.

# 'os.spawnl' Event

os.spawnl(os.P_NOWAIT, "/bin/echo", "echo", "Hello, World!")

# Classification: Unsafe: Spawning processes with untrusted or unsafe commands.

# 'os.posix_spawn' Event

os.posix_spawn("/bin/echo", ["/bin/echo", "Hello, World!"], os.environ)

# Classification: Unsafe: Spawning processes with untrusted or unsafe commands.

# 'os.kill' Event

import signal

os.kill(os.getpid(), signal.SIGTERM)

# Classification: Unsafe: Sending signal to terminate critical system processes or services.

# 'os.putenv' Event

os.putenv("MY_ENV_VAR", "value")

# Classification: Unsafe: Setting environment variables that could alter system behavior or expose sensitive data.

# Network Operations

# 'socket.connect' Event

import socket

s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s.connect(("www.example.com", 80))

# Classification: Unsafe: Connecting to unknown or untrusted servers, potentially exposing the system to malicious attacks.

# 'socket.sendto' Event

s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
s.sendto(b"Hello, World!", ("localhost", 9999))

# Classification: Unsafe: Sending data to untrusted or potentially malicious addresses.

# 'http.client.connect' Event

import http.client

conn = http.client.HTTPConnection("www.example.com")
conn.connect()

# Classification: Unsafe: Connecting to untrusted or potentially malicious servers.

# 'http.client.send' Event

conn = http.client.HTTPConnection("www.example.com")
conn.request("GET", "/")
response = conn.getresponse()
print(response.status, response.reason)

# Classification: Unsafe: Sending requests to untrusted or potentially malicious servers.

# 'ftplib.connect' Event

import ftplib

ftp = ftplib.FTP()
ftp.connect("ftp.example.com", 21)

# Classification: Unsafe: Connecting to untrusted FTP servers, potentially exposing credentials or data.

# 'smtplib.connect' Event

import smtplib

server = smtplib.SMTP()
server.connect("smtp.example.com", 25)

# Classification: Unsafe: Connecting to untrusted SMTP servers, potentially exposing email data.

# 'poplib.connect' Event

import poplib

pop_con = poplib.POP3("pop.example.com")

# Classification: Unsafe: Connecting to untrusted POP3 servers, potentially exposing email data.

# 'imaplib.open' Event

import imaplib

imap_conn = imaplib.IMAP4("imap.example.com")

# Classification: Unsafe: Connecting to untrusted IMAP servers, potentially exposing email data.

# 'urllib.Request' Event

import urllib.request

req = urllib.request.Request("http://www.example.com")
urllib.request.urlopen(req)

# Classification: Unsafe: Making requests to untrusted or potentially malicious web resources.

# Module Execution

# 'exec' Event

code = "print('Hello, World!')"
exec(code)

# Classification: Unsafe: Executing dynamically generated code from untrusted sources.

# 'compile' Event

code = "print('Hello, World!')"
compiled_code = compile(code, '<string>', 'exec')
exec(compiled_code)

# Classification: Unsafe: Compiling dynamically generated code from untrusted sources.

# 'code.__new__' Event

code_obj = compile("print('Hello, World!)", "<string>", "exec")

# Classification: Unsafe: Compiling dynamically generated code from untrusted sources.

# 'import' Event

import json

# Classification: Unsafe: Importing modules from untrusted sources or directories

# 'marshal.loads' Event

import marshal

data = marshal.dumps([1, 2, 3])
loaded_data = marshal.loads(data)

# Classification: Unsafe: Loading serialized data from untrusted sources, which could lead to arbitrary code execution.

# 'marshal.dumps' Event

data = marshal.dumps([1, 2, 3])
print(data)

# Classification: Unsafe: Serializing sensitive data that could be exposed.

# 'pickle.find_class' Event

import pickle

class MyClass:
    pass

data = pickle.dumps(MyClass())
loaded_obj = pickle.loads(data)

# Classification: Unsafe: Unpickling objects from untrusted sources, which could lead to arbitrary code execution.

# Object Management

# 'gc.get_objects' Event

import gc

gc.get_objects()

# Classification: Unsafe: Potential privacy or security issues if sensitive objects are exposed.

# 'rc.get_referents' Event

obj = [1, 2, 3]
gc.get_referents(obj)

# Classification: Unsafe: Potential privacy issues if sensitive referents are exposed.

# 'rc.get_referrers' Event

obj = [1, 2, 3]
print(gc.get_referrers(obj))

# Classification: Unsafe: Potential privacy issues if sensitive referrers are exposed.

# 'function.__new__' Event

def test_func():
    pass

# Classification: Unsafe: Creating functions that execute untrusted code

# 'object.__setattr__' Event

class TestClass:
    def __init__(self, name) -> None:
        self.name = name

obj = TestClass("Initial name")
obj.name = "Updated name"

# Classification: Unsafe: Setting attributes that may alter critical system behavior or expose sensitive data.

# 'object.__getattr__' Event

obj = TestClass("Example")
print(obj.name)

# Classification: Unsafe: Accessing attributes that may reveal sensitive data.

# 'object.__delattr__' Event

obj = TestClass("Example")
del obj.name

# Classification: Unsafe: Deleting attributes that may alter critical system behavior or expose sensitive data.

# 'ctypes.dlopen' Event

import ctypes

libc = ctypes.CDLL("libc.so.6")

# Classification: Unsafe: Loading untrusted or potentially malicious shared libraries.

# 'ctypes.set_errno' Event

ctypes.set_errno(0)

# Classification: Unsafe: Setting errno in ways that may alter critical system behavior or mask errors.
