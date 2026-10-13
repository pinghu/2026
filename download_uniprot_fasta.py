#!/usr/bin/env python3
"""Download FASTA sequences for UniRef90 identifiers in uniqUnirefID."""

import argparse
import concurrent.futures
import json
import sys
import time
from http.client import RemoteDisconnected
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen


BASE_URL = "https://rest.uniprot.org"
USER_AGENT = "uniqUnirefID-fasta-downloader/1.0"


def request_bytes(url, data=None, attempts=5):
    for attempt in range(attempts):
        request = Request(
            url,
            data=data,
            headers={
                "User-Agent": USER_AGENT,
                "Content-Type": "application/x-www-form-urlencoded" if data else "application/json",
            },
        )
        try:
            with urlopen(request, timeout=120) as response:
                return response.read()
        except HTTPError as error:
            if error.code == 404:
                return None
            if error.code not in (408, 429, 500, 502, 503, 504) or attempt + 1 == attempts:
                raise
            retry_after = error.headers.get("Retry-After")
            delay = float(retry_after) if retry_after and retry_after.isdigit() else 2**attempt
            time.sleep(delay)
        except (TimeoutError, URLError, RemoteDisconnected, ConnectionResetError):
            if attempt + 1 == attempts:
                raise
            time.sleep(2**attempt)
    return None


def read_ids(path):
    ids = []
    seen = set()
    with path.open(encoding="utf-8") as source:
        for line_number, line in enumerate(source, start=1):
            value = line.strip()
            if not value or value.startswith("#") or value == "UNMAPPED":
                continue
            if not value.startswith("UniRef90_"):
                raise ValueError(f"{path}:{line_number}: expected a UniRef90_ identifier, got {value!r}")
            if value not in seen:
                ids.append(value)
                seen.add(value)
    return ids


def map_to_uniparc(accessions, poll_seconds):
    body = urlencode(
        {
            "from": "UniProtKB_AC-ID",
            "to": "UniParc",
            "ids": ",".join(accessions),
        }
    ).encode()
    response = request_bytes(
        f"{BASE_URL}/idmapping/run",
        data=body,
    )
    job_id = json.loads(response)["jobId"]

    while True:
        status_bytes = request_bytes(f"{BASE_URL}/idmapping/status/{job_id}")
        status = json.loads(status_bytes)
        if "results" in status:
            results_bytes = request_bytes(f"{BASE_URL}/idmapping/stream/{job_id}?format=json")
            rows = json.loads(results_bytes)["results"]
            break
        job_status = status.get("jobStatus")
        if job_status in ("NEW", "RUNNING"):
            time.sleep(poll_seconds)
            continue
        if job_status == "FINISHED":
            results_bytes = request_bytes(f"{BASE_URL}/idmapping/stream/{job_id}?format=json")
            rows = json.loads(results_bytes)["results"]
            break
        raise RuntimeError(f"UniProt ID mapping job {job_id} ended unexpectedly: {status}")

    mapped = {}
    for row in rows:
        source = row["from"]
        target = row["to"]
        if isinstance(target, dict):
            target = target.get("uniParcId")
        if source and target:
            mapped[source] = target
    return mapped


def fetch_uniparc_fasta(uniparc_ids):
    query = " OR ".join(f"upi:{uniparc_id}" for uniparc_id in uniparc_ids)
    params = urlencode({"query": query, "format": "fasta", "size": len(uniparc_ids)})
    fasta_bytes = request_bytes(f"{BASE_URL}/uniparc/search?{params}")
    sequences = {}
    if fasta_bytes is None:
        return sequences

    current_id = None
    sequence_parts = []
    for line in fasta_bytes.decode("utf-8").splitlines():
        if line.startswith(">"):
            if current_id and sequence_parts:
                sequences[current_id] = "".join(sequence_parts)
            current_id = line[1:].split()[0]
            sequence_parts = []
        elif current_id:
            sequence_parts.append(line.strip())
    if current_id and sequence_parts:
        sequences[current_id] = "".join(sequence_parts)
    return sequences


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("uniqUnirefID"))
    parser.add_argument("--output", type=Path, default=Path("uniprot_sequences.fasta"))
    parser.add_argument("--missing", type=Path, default=Path("uniprot_missing_ids.txt"))
    parser.add_argument("--batch-size", type=int, default=10000)
    parser.add_argument("--fasta-batch-size", type=int, default=40)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--poll-seconds", type=int, default=5)
    parser.add_argument("--limit", type=int, help="Process only the first N unique identifiers")
    parser.add_argument("--resume", action="store_true", help="Append to existing outputs and skip completed IDs")
    args = parser.parse_args()
    if min(args.batch_size, args.fasta_batch_size, args.workers, args.poll_seconds) < 1:
        parser.error("batch sizes, workers, and poll interval must be positive")

    ids = read_ids(args.input)
    if args.limit is not None:
        if args.limit < 1:
            parser.error("--limit must be at least 1")
        ids = ids[: args.limit]
    original_id_by_accession = {
        uniref_id.removeprefix("UniRef90_"): uniref_id for uniref_id in ids
    }
    found = set()
    missing = set()
    if args.resume:
        if args.output.exists():
            with args.output.open(encoding="utf-8") as fasta_in:
                found = {
                    line[1:].split()[0]
                    for line in fasta_in
                    if line.startswith(">")
                }
        if args.missing.exists():
            with args.missing.open(encoding="utf-8") as missing_in:
                missing = {line.strip() for line in missing_in if line.strip()}

    output_mode = "a" if args.resume else "w"
    with args.output.open(output_mode, encoding="utf-8") as fasta_out:
        for start in range(0, len(ids), args.batch_size):
            batch_ids = [
                uniref_id
                for uniref_id in ids[start : start + args.batch_size]
                if uniref_id not in found and uniref_id not in missing
            ]
            if not batch_ids:
                continue
            accessions = [uniref_id.removeprefix("UniRef90_") for uniref_id in batch_ids]
            mapped = map_to_uniparc(accessions, args.poll_seconds)
            mapped_ids = [(accession, mapped[accession]) for accession in accessions if accession in mapped]
            missing.update(
                original_id_by_accession[accession]
                for accession in accessions
                if accession not in mapped
            )
            sequence_batches = [
                mapped_ids[index : index + args.fasta_batch_size]
                for index in range(0, len(mapped_ids), args.fasta_batch_size)
            ]
            with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
                futures = {
                    pool.submit(fetch_uniparc_fasta, [upi for _, upi in sequence_batch]): sequence_batch
                    for sequence_batch in sequence_batches
                }
                for future in concurrent.futures.as_completed(futures):
                    sequence_batch = futures[future]
                    sequences = future.result()
                    for accession, upi in sequence_batch:
                        sequence = sequences.get(upi)
                        if not sequence:
                            missing.add(original_id_by_accession[accession])
                            continue
                        uniref_id = original_id_by_accession[accession]
                        fasta_out.write(f">{uniref_id} source=UniParc:{upi}\n")
                        for offset in range(0, len(sequence), 60):
                            fasta_out.write(sequence[offset : offset + 60] + "\n")
                        found.add(uniref_id)

            with args.missing.open("w", encoding="utf-8") as missing_out:
                missing_out.writelines(f"{uniref_id}\n" for uniref_id in sorted(missing))
            fasta_out.flush()
            print(
                f"Processed {len(found) + len(missing)}/{len(ids)} "
                f"(sequences: {len(found)}, missing: {len(missing)})",
                file=sys.stderr,
                flush=True,
            )

    print(f"Processed {len(ids)} unique identifiers; wrote {len(found)} sequences and {len(missing)} missing IDs.")
    print(f"FASTA: {args.output}")
    print(f"Missing IDs: {args.missing}")


if __name__ == "__main__":
    main()
