{
  description = "Generic Python project dev shell with venv and requirements.txt support";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";

  outputs = { self, nixpkgs }:
    let
      system = "x86_64-linux";  # or "aarch64-linux" for ARM
      pkgs = import nixpkgs { inherit system; };
    in {
      devShells.${system}.default = pkgs.mkShell {
        buildInputs = [
          pkgs.python3  # change to pkgs.python311, pkgs.python312, etc.
          pkgs.python3Packages.virtualenv
        ];

        shellHook = ''
          # Create or activate venv
          if [ ! -d .venv ]; then
            echo "🔧 Creating new virtual environment..."
            python -m venv .venv
            source .venv/bin/activate
            pip install --upgrade pip
            if [ -f requirements.txt ]; then
              echo "📦 Installing requirements.txt..."
              pip install -r requirements.txt
            fi
          else
            source .venv/bin/activate
          fi

          echo "🐍 Python $(python --version) with $(pip list | wc -l) packages ready."
          echo "Run: pytest"
        '';
      };
    };
}

