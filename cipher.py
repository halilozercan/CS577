import os
import hashlib
import Crypto.Cipher.AES as AES


class AesCtrCipher:
    @staticmethod
    def md5sum(raw):
        m = hashlib.md5()
        m.update(raw)
        return m.hexdigest()

    @staticmethod
    def pad(s):
        """note that the padding is no necessary"""
        """return s + (Cipher.BS - len(s) % Cipher.BS) * chr(Cipher.BS - len(s) % Cipher.BS)"""
        return s

    @staticmethod
    def unpad(s):
        """return s[0:-ord(s[-1])]"""
        return s

    def __init__(self, key):
        self.key = AesCtrCipher.md5sum(key)
        # the state of the counter callback
        self.cnter_cb_called = 0
        self.secret = None

    def _reset_counter_callback_state(self, secret):
        self.cnter_cb_called = 0
        self.secret = secret

    def _counter_callback(self):
        """
                this function should be stateful
                """
        self.cnter_cb_called += 1
        return self.secret[self.cnter_cb_called % AES.block_size] * AES.block_size

    def encrypt(self, raw):
        secret = os.urandom(AES.block_size)  # random choose a "secret" which is not secret
        self._reset_counter_callback_state(secret)
        cipher = AES.new(self.key, AES.MODE_CTR, counter=self._counter_callback)
        raw_padded = AesCtrCipher.pad(raw)
        enc_padded = cipher.encrypt(raw_padded)
        return secret + enc_padded  # yes, it is not secret

    def decrypt(self, enc):
        secret = enc[:AES.block_size]
        self._reset_counter_callback_state(secret)
        cipher = AES.new(self.key, AES.MODE_CTR, counter=self._counter_callback)
        enc_padded = enc[AES.block_size:]  # we didn't encrypt the secret, so don't decrypt it
        raw_padded = cipher.decrypt(enc_padded)
        return AesCtrCipher.unpad(raw_padded)
