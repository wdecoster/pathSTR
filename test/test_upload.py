import base64
import gzip
import unittest

from pathSTR import parse_input
from pathSTR.plot import hex_to_rgb


def data_url(payload: bytes) -> str:
    return "data:application/gzip;base64," + base64.b64encode(payload).decode()


class DecodeUploadTests(unittest.TestCase):
    def test_gzipped_vcf_is_decompressed(self):
        with open("test/HG00641.vcf.gz", "rb") as f:
            raw = f.read()
        decoded = parse_input.decode_upload(data_url(raw))
        self.assertTrue(decoded.startswith(b"##fileformat=VCF"))

    def test_plain_content_is_passed_through(self):
        self.assertEqual(parse_input.decode_upload(data_url(b"##fileformat")), b"##fileformat")

    def test_malformed_content_is_rejected(self):
        for contents in [None, "no-comma", "data:x;base64,@@not base64@@"]:
            with self.assertRaises(parse_input.UploadError):
                parse_input.decode_upload(contents)

    def test_oversized_upload_is_rejected(self):
        with self.assertRaises(parse_input.UploadError):
            parse_input.decode_upload(data_url(b"A" * (parse_input.MAX_UPLOAD_BYTES + 10)))

    def test_gzip_bomb_is_rejected(self):
        bomb = gzip.compress(b"\0" * (parse_input.MAX_DECOMPRESSED_BYTES + 1000))
        self.assertLess(len(bomb), parse_input.MAX_UPLOAD_BYTES)
        with self.assertRaises(parse_input.UploadError):
            parse_input.decode_upload(data_url(bomb))

    def test_parse_uploaded_vcf_returns_none_on_garbage(self):
        self.assertIsNone(
            parse_input.parse_uploaded_vcf("garbage", "x.vcf.gz", None, "hg38", "STRdust")
        )


class ColorTests(unittest.TestCase):
    def test_hex_to_rgb(self):
        self.assertEqual(hex_to_rgb("#1f77b4"), "rgb(31, 119, 180)")


if __name__ == "__main__":
    unittest.main()
