import os
import datetime
import gurobipy as gp

class Config:
    def __init__(self, output_dir=None, output_name=None, save_lp=False, save_sol=False, save_ilp=False):
        """
        Configuration for output file management.

        :param output_dir:  User-specified output directory (defaults to the current working directory)
        :param output_name: User-specified file name prefix (defaults to 'YYYYMMDD_HHMM' format)
        :param save_lp:     Whether to save the LP file
        :param save_sol:    Whether to save the SOL file
        :param save_ilp:    Whether to save the ILP file
        """
        # If the output directory is not specified, use the current working directory
        self.output_dir = output_dir if output_dir else os.getcwd()
        
        # If the output name is not specified, use the current timestamp in 'YYYYMMDD_HHMM' format
        if output_name:
            self.output_name = output_name
        else:
            now = datetime.datetime.now()
            self.output_name = now.strftime("%Y%m%d_%H%M")
        
        # File output settings
        self.save_lp = save_lp
        self.save_sol = save_sol
        self.save_ilp = save_ilp

        # Ensure the output directory exists
        os.makedirs(self.output_dir, exist_ok=True)

    def get_output_path(self, file_type):
        """Returns the full output file path for a given file type (LP, SOL, or ILP)."""
        return os.path.join(self.output_dir, f"{self.output_name}.{file_type}")

def save_model_info(m: gp.Model, cfg: Config):
    """
    Saves model files based on the configuration.

    :param m:   Gurobi model instance
    :param cfg: Config instance specifying which files to save
    """
    if cfg.save_lp:
        m.write(cfg.get_output_path("lp"))
    if cfg.save_sol:
        m.write(cfg.get_output_path("sol"))
    if cfg.save_ilp:
        m.computeIIS()
        m.write(cfg.get_output_path("ilp"))
