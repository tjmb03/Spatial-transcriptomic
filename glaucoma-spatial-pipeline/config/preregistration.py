# preregistration.py — placeholder until hash is locked before data collection
PREREG_HASH = "UNSET"

PREREG = {
    "version": "1.0.0",
    "registration_doi": "UNSET",
    "registration_date": "UNSET",
}

def compute_hash():
    import hashlib, json
    return hashlib.sha256(
        json.dumps(PREREG, sort_keys=True).encode()
    ).hexdigest()

def validate_preregistration(strict=False):
    if PREREG_HASH == "UNSET":
        print("[PREREG] WARNING: Hash not set — run in development mode only")
        return
    current = compute_hash()
    if current != PREREG_HASH:
        msg = f"[PREREG] HASH MISMATCH\n  Expected: {PREREG_HASH}\n  Current:  {current}"
        if strict:
            raise ValueError(msg)
        print(msg)
    else:
        print("[PREREG] Hash verified OK")
