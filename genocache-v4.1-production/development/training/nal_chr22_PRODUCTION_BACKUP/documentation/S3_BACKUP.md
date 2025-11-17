# S3 Backup Instructions

This guide shows how to backup the Chr22 NAL training files to S3.

## Prerequisites

1. Install AWS CLI (if not already installed):
```bash
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip
sudo ./aws/install
```

2. Configure AWS credentials:
```bash
aws configure
# Enter:
# - AWS Access Key ID
# - AWS Secret Access Key  
# - Default region (e.g., us-east-1)
# - Default output format (json)
```

## Quick Backup Commands

### 1. Backup Entire Folder to S3

```bash
cd /home/nebius/genocache/genocache-v4.1-production/development/training

# Sync entire backup folder (excludes models/ due to size)
aws s3 sync nal_chr22_PRODUCTION_BACKUP/ \
  s3://YOUR-BUCKET-NAME/nal_chr22_backup/ \
  --exclude "models/*pt" \
  --exclude "indexes/*"
```

### 2. Backup Specific Files Only

```bash
# Backup training code
aws s3 cp nal_chr22_PRODUCTION_BACKUP/training/ \
  s3://YOUR-BUCKET-NAME/nal_chr22_backup/training/ \
  --recursive

# Backup alignment code
aws s3 cp nal_chr22_PRODUCTION_BACKUP/alignment/ \
  s3://YOUR-BUCKET-NAME/nal_chr22_backup/alignment/ \
  --recursive

# Backup documentation
aws s3 cp nal_chr22_PRODUCTION_BACKUP/documentation/ \
  s3://YOUR-BUCKET-NAME/nal_chr22_backup/documentation/ \
  --recursive
```

### 3. Backup Models Separately (Optional)

Models are large (5.5 MB each), so backup separately:

```bash
# If models exist locally
aws s3 cp ../nal_single_chr/models/chr22_nal_512bp_final.pt \
  s3://YOUR-BUCKET-NAME/nal_chr22_backup/models/
```

### 4. Backup Indexes (Optional)

Indexes are LARGE (1.1 GB), backup only if needed:

```bash
aws s3 cp ../nal_single_chr/indexes/ \
  s3://YOUR-BUCKET-NAME/nal_chr22_backup/indexes/ \
  --recursive
```

## Restore From S3

### Restore Code Only (Fast)

```bash
aws s3 sync s3://YOUR-BUCKET-NAME/nal_chr22_backup/ \
  ./nal_chr22_RESTORED/ \
  --exclude "models/*" \
  --exclude "indexes/*"
```

### Restore Everything (Including Models & Indexes)

```bash
aws s3 sync s3://YOUR-BUCKET-NAME/nal_chr22_backup/ \
  ./nal_chr22_RESTORED/
```

## What to Backup

### Essential (Small, ~1 MB):
- ✅ `training/*.py` - Training scripts
- ✅ `indexing/*.py` - Index building
- ✅ `alignment/*.py` - Alignment pipeline
- ✅ `documentation/*.md` - All docs
- ✅ `results/*.log` - Test results
- ✅ `README.md` - Main readme

### Optional (Large):
- ⚠️  `models/*.pt` - 5.5 MB each (can retrain)
- ⚠️  `indexes/*` - 1.1 GB total (can rebuild)

## Size Summary

```
Code files:        ~100 KB
Documentation:     ~50 KB
Results:           ~10 KB
Models:            ~5.5 MB (optional)
Indexes:           ~1.1 GB (optional)
```

## Recommended Backup Strategy

**For code backup (recommended):**
```bash
# Backup just code and docs (~150 KB)
aws s3 sync nal_chr22_PRODUCTION_BACKUP/ \
  s3://YOUR-BUCKET-NAME/nal_chr22_backup/ \
  --exclude "models/*" \
  --exclude "indexes/*" \
  --exclude "*.pyc" \
  --exclude "__pycache__/*"
```

**For complete backup (if needed):**
```bash
# Backup everything including models
aws s3 sync nal_chr22_PRODUCTION_BACKUP/ \
  s3://YOUR-BUCKET-NAME/nal_chr22_backup/

# Also backup models from original location
aws s3 sync ../nal_single_chr/models/ \
  s3://YOUR-BUCKET-NAME/nal_chr22_backup/models/

# And indexes
aws s3 sync ../nal_single_chr/indexes/ \
  s3://YOUR-BUCKET-NAME/nal_chr22_backup/indexes/
```

## Verification

Check what was uploaded:
```bash
aws s3 ls s3://YOUR-BUCKET-NAME/nal_chr22_backup/ --recursive --human-readable
```

## Cost Estimate

- **Code only** (~150 KB): $0.000003/month
- **Code + Models** (~6 MB): $0.00014/month  
- **Everything** (~1.1 GB): $0.025/month

(Based on S3 Standard storage pricing)

## Notes

- Replace `YOUR-BUCKET-NAME` with your actual S3 bucket name
- Consider using `--dryrun` flag first to see what will be copied
- Use `--delete` flag with `sync` to remove files from S3 that were deleted locally
- For large files, consider using `--storage-class GLACIER` for cheaper long-term storage
