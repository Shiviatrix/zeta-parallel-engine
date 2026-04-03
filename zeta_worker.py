import requests
import subprocess
import time
import os

SERVER_URL = "http://127.0.0.1:5001" # Change to Cloud IP if remote
ZETA_BIN = "./zeta_v7"

def run_worker():
    while True:
        try:
            # 1. Request Work
            print("Requesting work...")
            resp = requests.get(f"{SERVER_URL}/get_work")
            if resp.status_code != 200:
                print("Server error. Retrying in 5s...")
                time.sleep(5)
                continue
            
            job = resp.json()
            job_id = job["job_id"]
            start_n = job["start_n"]
            size = job["size"]
            
            print(f"Received Job {job_id}: Start={start_n}, Size={size}")
            
            # 2. Run Engine (Batch Mode)
            # ./zeta_v6 <start_n> <max_zeros>
            cmd = [ZETA_BIN, str(start_n), str(size)]
            
            # Use subprocess to capture output
            start_time = time.time()
            process = subprocess.run(cmd, capture_output=True, text=True)
            duration = time.time() - start_time
            
            # Parse Output to confirm verified count (simplified)
            # Real implementation would parse "Verified: X" from stdout if needed
            # For now, assume if exit code 0, all requested zeroes were verified (batch limit reached)
            verified_count = size 
            
            if process.returncode != 0:
                print(f"Engine Failed: {process.stderr}")
                verified_count = 0
            
            print(f"Job Finished in {duration:.2f}s.")
            
            # 3. Submit Results
            payload = {
                "job_id": job_id,
                "verified_count": verified_count,
                "duration": duration
            }
            requests.post(f"{SERVER_URL}/submit_work", json=payload)
            print(f"Results submitted. (Speed: {verified_count/duration:.2f} z/s)\n")
            
        except Exception as e:
            print(f"Worker Error: {e}")
            time.sleep(5)

if __name__ == "__main__":
    if not os.path.exists(ZETA_BIN):
        print(f"Error: {ZETA_BIN} not found. Compile it first!")
    else:
        run_worker()
