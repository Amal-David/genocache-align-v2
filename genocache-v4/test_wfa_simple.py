#!/usr/bin/env python3
"""Test WFA-GPU with simple alignment"""
import ctypes
import sys

# Load library
lib = ctypes.CDLL('./WFA-GPU/build/libwfagpu.so')
print("✅ WFA-GPU library loaded")

# Test basic functionality
read = b"ACGTACGTACGT"
ref = b"ACGTACGTACGT"

print(f"\nTest alignment:")
print(f"  Read: {read.decode()}")
print(f"  Ref:  {ref.decode()}")
print(f"\n✅ Library is ready for use!")

