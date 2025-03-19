import importlib


def select_rocket(rocket_name: str):
    try:
        # Dynamically import the rocket module
        rocket_module = importlib.import_module(f'launch_vehicles.{rocket_name}')
        return rocket_module
    except ModuleNotFoundError:
        raise ValueError(f'Rocket {rocket_name} not found in launch_vehicles.')