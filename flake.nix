{
  description = "A Nix-flake-based and Snakemake based bioinformatics pipeline development environment";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-24.11";
    flake-utils.url = "github:numtide/flake-utils";
    git-hooks.url = "github:cachix/git-hooks.nix";
  };

  outputs = {
    self,
    nixpkgs,
    flake-utils,
    git-hooks,
  }:
    flake-utils.lib.eachDefaultSystem (system: let
      overlays = [
        (self: super: rec {
          python = super.python;
          snakemake = super.pythonPackages.snakemake;
        })
      ];
      pkgs = import nixpkgs {inherit overlays system;};
      packages = with pkgs; [
        prettier

        git
        typos
        alejandra
      ];
    in {
      checks = {
        pre-commit-check = git-hooks.lib.${system}.run {
          src = ./.;
          hooks = {
            typos.enable = true; # Source code spell checker
            alejandra.enable = true; # Nix linter
            prettier.enable = true; # Markdown & TS formatter
            nixpkgs-fmt.enable = true;

            snakefmt = {
              enable = true;
              name = "snakefmt";
              description = "Run snakefmt";
              files = "\.smk$";
              entry = "${pkgs.snakefmt}/bin/snakefmt";
            };
          };
          settings = {
            typos = {
              write = true; # Automatically fix typos
              ignored-words = [];
            };
            prettier = {
              write = true; # Automatically format files
              configPath = "./.prettierrc.yaml";
            };
          };
        };
      };

      devShells.default = pkgs.mkShell {
        inherit packages;

        shellHook = ''
          ${self.checks.${system}.pre-commit-check.shellHook}
        '';
      };
    });
}
