import unittest
from unittest.mock import Mock, patch

import app


class AppEntrypointTests(unittest.TestCase):
    def test_wsgi_factory_uses_database_environment(self):
        dash_app = Mock()
        with patch.dict("os.environ", {"PATHSTR_DB": "/mounted/database"}):
            with patch.object(app, "create_app", return_value=dash_app) as create:
                server = app.create_wsgi_app()

        self.assertIs(server, dash_app.server)
        args = create.call_args.args[0]
        self.assertEqual(args.db, "/mounted/database")
        self.assertFalse(args.store_only)
        self.assertEqual(create.call_args.kwargs, {"log_to_file": False})

    def test_wsgi_factory_default_database(self):
        with patch.dict("os.environ", {}, clear=True):
            with patch.object(app, "create_app", return_value=Mock()) as create:
                app.create_wsgi_app()
        self.assertEqual(create.call_args.args[0].db, "/data/pathSTR_database")

    def test_cli_starts_dash_with_requested_options(self):
        args = app.get_args(["--db", "db", "--host", "127.0.0.1", "--port", "9000", "--debug"])
        dash_app = Mock()
        with patch.object(app, "get_args", return_value=args):
            with patch.object(app, "create_app", return_value=dash_app):
                app.main()
        dash_app.run.assert_called_once_with(host="127.0.0.1", port=9000, debug=True)

    def test_cli_refuses_debug_on_network_host(self):
        args = app.get_args(["--db", "db", "--host", "0.0.0.0", "--debug"])
        with patch.object(app, "get_args", return_value=args):
            with patch.object(app, "create_app") as create:
                with self.assertRaises(SystemExit):
                    app.main()
        create.assert_not_called()

    def test_cli_store_only_does_not_start_server(self):
        args = app.get_args(["--db", "db", "--store_only"])
        with patch.object(app, "get_args", return_value=args):
            with patch.object(app, "create_app", return_value=None) as create:
                app.main()
        create.assert_called_once_with(args)


if __name__ == "__main__":
    unittest.main()
