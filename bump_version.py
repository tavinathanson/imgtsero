#!/usr/bin/env python3
"""Version bumping script for imgtsero."""

import argparse
import re
import subprocess
import sys
from pathlib import Path


def update_file(filepath, pattern, replacement):
    """Update version in a file."""
    content = filepath.read_text()
    new_content = re.sub(pattern, replacement, content)
    if content == new_content:
        print(f"Warning: No changes made to {filepath}")
        return False
    filepath.write_text(new_content)
    return True


def main():
    parser = argparse.ArgumentParser(description="Bump imgtsero version")
    parser.add_argument("version", help="New version number (e.g., 0.4.1)")
    parser.add_argument("-m", "--message", help="Git commit and tag message", 
                        default="Release v{version}")
    parser.add_argument("--no-git", action="store_true", 
                        help="Skip git operations (commit and tag)")
    args = parser.parse_args()

    # Validate version format
    if not re.match(r'^\d+\.\d+\.\d+$', args.version):
        print(f"Error: Invalid version format '{args.version}'. Use X.Y.Z format.")
        sys.exit(1)

    # Update files
    root = Path(__file__).parent
    files_updated = 0
    
    # Update __init__.py
    init_file = root / "imgtsero" / "__init__.py"
    if update_file(init_file, 
                   r'__version__ = "[^"]*"', 
                   f'__version__ = "{args.version}"'):
        print(f"✓ Updated {init_file}")
        files_updated += 1

    # Update pyproject.toml
    pyproject_file = root / "pyproject.toml"
    if update_file(pyproject_file, 
                   r'version = "[^"]*"', 
                   f'version = "{args.version}"'):
        print(f"✓ Updated {pyproject_file}")
        files_updated += 1

    # Update setup.py
    setup_file = root / "setup.py"
    if update_file(setup_file, 
                   r'version="[^"]*"', 
                   f'version="{args.version}"'):
        print(f"✓ Updated {setup_file}")
        files_updated += 1

    if files_updated == 0:
        print("No files were updated. Version may already be set.")
        sys.exit(1)

    print(f"\nVersion bumped to {args.version} in {files_updated} files")

    # Git operations
    if not args.no_git:
        try:
            # Stage the changed files
            subprocess.run(["git", "add", "imgtsero/__init__.py", "pyproject.toml", "setup.py"], 
                         check=True)
            
            # Commit
            commit_msg = f"Bump version to {args.version}"
            subprocess.run(["git", "commit", "-m", commit_msg], check=True)
            print(f"✓ Committed: {commit_msg}")
            
            # Create tag
            tag_name = f"v{args.version}"
            tag_msg = args.message.format(version=args.version)
            subprocess.run(["git", "tag", "-a", tag_name, "-m", tag_msg], check=True)
            print(f"✓ Created tag: {tag_name}")
            
            print(f"\nNext steps:")
            print(f"  git push origin main")
            print(f"  git push origin {tag_name}")
            
        except subprocess.CalledProcessError as e:
            print(f"\nError during git operations: {e}")
            print("Files have been updated but not committed.")
            sys.exit(1)


if __name__ == "__main__":
    main()