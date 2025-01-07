import pytest #type:ignore
from dash.testing.application_runners import import_app #type:ignore
import yaml #type:ignore

class TestDashApp:
    
    """ Setup instance of dash app for testing. """
    @pytest.fixture(scope="class")
    def setup(self, dash_duo):
        self.dash_duo = dash_duo
        self.dash_duo.start_server(import_app("src.dash_app"))
        self.config_path = "config.yaml"
        
    def test_save_config_button(self):
        population_slider = self.dash_duo.find_element("#population-size-slider")
        # Clear and set population slider to 100
        self.dash_duo.clear_input(population_slider)
        population_slider.send_keys("100")
        
        # Click the Save Config button
        save_config_button = self.dash_duo.find_element("#save-button")
        save_config_button.click()
        
        # TODO - assert that population has been set properly in config
        with open(self.config_path, "r") as file:
            current_config = yaml.safe_load(file)
        assert current_config["population_size"] == 200
    
    def test_start_ga_button(self):
        # Click the Start GA button
        start_button = self.dash_duo.find_element("#start-ga-buton")
        start_button.click()
        
        ga_status = self.dash_duo.find_element("#ga-status").text
        assert ga_status == "GA Started..."
        
    def test_update_results_button(self):
        # Click the Update Results button
        update_results_button = self.dash_duo.find_element("#update-results-button")
        update_results_button.click()
        
        # Verify diversity graph is loaded
        diversity_graph = self.dash_duo.find_element("#diversity-graph")
        assert diversity_graph is not None
        
        # TODO - extend to other graphs/results