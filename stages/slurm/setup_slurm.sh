#!/bin/bash
# stages/slurm/setup_slurm.sh — Install and configure single-node SLURM on Debian
# Sets up a minimal SLURM environment for benchmarking on the local machine.
#
# Usage: sudo bash stages/slurm/setup_slurm.sh
#
# This creates a single-node SLURM cluster suitable for scaling experiments.
# After running, verify with: sinfo -N -l

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Detect CPU and memory
NCPUS=$(nproc)
TOTAL_MEM_MB=$(free -m | awk '/^Mem:/{print $2}')
HOSTNAME=$(hostname -s)

echo "=== SLURM Single-Node Setup ==="
echo "  Host:   $HOSTNAME"
echo "  CPUs:   $NCPUS"
echo "  Memory: ${TOTAL_MEM_MB}MB"
echo ""

# 1. Install packages
echo "[1/5] Installing SLURM packages..."
apt-get update -qq
apt-get install -y -qq slurm-wlm slurm-client munge libmunge-dev

# 2. Configure Munge (authentication)
echo "[2/5] Configuring Munge..."
if [[ ! -f /etc/munge/munge.key ]]; then
    dd if=/dev/urandom bs=1 count=1024 > /etc/munge/munge.key 2>/dev/null
    chown munge:munge /etc/munge/munge.key
    chmod 400 /etc/munge/munge.key
fi
systemctl enable munge
systemctl restart munge

# 3. Generate slurm.conf
echo "[3/5] Generating slurm.conf..."
cat > /etc/slurm/slurm.conf << SLURM_CONF
# SLURM configuration for single-node benchmarking
ClusterName=edna-benchmark
SlurmctldHost=${HOSTNAME}

# Scheduler
SchedulerType=sched/backfill
SelectType=select/cons_tres
SelectTypeParameters=CR_Core_Memory

# Logging
SlurmctldLogFile=/var/log/slurm/slurmctld.log
SlurmdLogFile=/var/log/slurm/slurmd.log
SlurmctldPidFile=/run/slurmctld.pid
SlurmdPidFile=/run/slurmd.pid

# Process tracking (cgroup preferred on modern systems; linuxproc as fallback)
ProctrackType=proctrack/cgroup

# Timeouts
SlurmctldTimeout=120
SlurmdTimeout=120
InactiveLimit=0
MinJobAge=300
MaxJobCount=10000
ReturnToService=1

# Accounting (minimal — needed for job tracking)
AccountingStorageType=accounting_storage/none

# Partitions and nodes
NodeName=${HOSTNAME} CPUs=${NCPUS} RealMemory=$((TOTAL_MEM_MB - 512)) State=UNKNOWN
PartitionName=main Nodes=${HOSTNAME} Default=YES MaxTime=24:00:00 State=UP
SLURM_CONF

# 4. Create required directories
echo "[4/5] Creating directories and setting permissions..."
mkdir -p /var/log/slurm /var/spool/slurm/ctld /var/spool/slurm/d
chown -R slurm:slurm /var/log/slurm /var/spool/slurm 2>/dev/null || true

# 5. Start services
echo "[5/5] Starting SLURM services..."
systemctl enable slurmctld slurmd
systemctl restart slurmctld
sleep 2
systemctl restart slurmd
sleep 2

# Verify
echo ""
echo "=== Verification ==="
sinfo -N -l 2>&1 || echo "Warning: sinfo failed — SLURM may need a moment to start"
echo ""
scontrol show partition main 2>&1 || true
echo ""
echo "SLURM single-node setup complete."
echo "Run 'sinfo' to verify, then use 'sbatch' to submit jobs."
