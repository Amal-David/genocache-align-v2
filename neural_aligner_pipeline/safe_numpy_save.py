#!/usr/bin/env python3
"""
Safely save large numpy arrays with progress indication
"""
import os, sys
import numpy as np

def safe_save_large_array(array, filename, chunk_size=100_000):
    """
    Save large numpy array safely with chunking and verification
    """
    print(f"Saving {filename}...")
    print(f"  Shape: {array.shape}")
    print(f"  Size: {array.nbytes / 1024**3:.2f} GB")
    
    # Method 1: Try direct save with memmap
    temp_file = filename + ".tmp"
    
    try:
        print(f"  Using memory-mapped save...")
        # Create memory-mapped array
        mmap_array = np.lib.format.open_memmap(
            temp_file, 
            mode='w+', 
            dtype=array.dtype, 
            shape=array.shape
        )
        
        # Copy in chunks
        n_chunks = (len(array) + chunk_size - 1) // chunk_size
        for i in range(n_chunks):
            start_idx = i * chunk_size
            end_idx = min((i + 1) * chunk_size, len(array))
            mmap_array[start_idx:end_idx] = array[start_idx:end_idx]
            
            if (i + 1) % 10 == 0 or i == n_chunks - 1:
                progress = (end_idx / len(array)) * 100
                print(f"    Progress: {end_idx:,}/{len(array):,} ({progress:.1f}%)")
        
        # Flush to disk
        del mmap_array
        
        # Rename to final
        if os.path.exists(filename):
            os.remove(filename)
        os.rename(temp_file, filename)
        
        print(f"  ✓ Saved successfully")
        return True
        
    except Exception as e:
        print(f"  ✗ Memmap save failed: {e}")
        
        # Cleanup temp file
        if os.path.exists(temp_file):
            os.remove(temp_file)
        
        # Method 2: Fallback to standard numpy save
        try:
            print(f"  Trying standard numpy save...")
            np.save(filename, array)
            print(f"  ✓ Saved successfully")
            return True
        except Exception as e2:
            print(f"  ✗ Standard save also failed: {e2}")
            return False

def verify_saved_file(filename, expected_shape):
    """Verify saved numpy file"""
    print(f"Verifying {filename}...")
    try:
        # Load with memmap (doesn't load into RAM)
        data = np.load(filename, mmap_mode='r')
        
        if data.shape != expected_shape:
            print(f"  ✗ Shape mismatch: {data.shape} vs {expected_shape}")
            return False
        
        print(f"  ✓ Shape correct: {data.shape}")
        print(f"  ✓ File size: {os.path.getsize(filename) / 1024**3:.2f} GB")
        return True
        
    except Exception as e:
        print(f"  ✗ Verification failed: {e}")
        return False

if __name__ == "__main__":
    # Test with dummy data
    print("Testing safe save with dummy array...")
    test_array = np.random.randn(1000, 256).astype('float32')
    
    if safe_save_large_array(test_array, "test_save.npy"):
        if verify_saved_file("test_save.npy", test_array.shape):
            print("\n✓ Safe save system working!")
            os.remove("test_save.npy")
        else:
            print("\n✗ Verification failed")
    else:
        print("\n✗ Save failed")
