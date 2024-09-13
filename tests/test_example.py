import unittest
from unittest.mock  import patch
from automated_fits import main

class TestAutomatedFits(unittest.TestCase):

    # To test individual functions (if possible)
    #def test_single_function(self):
    #    result = single_function(5, 10)
    #    self.assertEqual(result, 50)  # Change this by the expected value

    @patch('sys.argv', ['automated_fits.py', '--x', '5', '--y', '10'])
    def test_main(self):
        # Prueba de la función main, simulando argumentos de línea de comandos
        resultado = main()
        self.assertEqual(resultado, 50)  
    def test_main(self):
        # Verify that ends without errors
        try:
            main()
        except Exception as e:
            self.fail(f"main() threw an exception: {e}")

        # Validate some kind of result?

if __name__ == '__main__':
    unittest.main()
